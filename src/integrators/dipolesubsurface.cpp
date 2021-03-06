
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// integrators/dipolesubsurface.cpp*
#include "stdafx.h"
#include "integrators/dipolesubsurface.h"
#include "scene.h"
#include "montecarlo.h"
#include "sampler.h"
#include "progressreporter.h"
#include "intersection.h"
#include "paramset.h"
#include "reflection.h"
#include "octree.h"
#include "camera.h"
#include "floatfile.h"

#define NUM_DIPOLE 21 //should be odd!
#define MAX_RADIUS 20.0f
#define D12_SIZE 2000
#define RADIUS_PREC 0.01f
#define THICKNESS 0.25
struct DiffusionReflectance;
struct MultipoleDR;

// DipoleSubsurfaceIntegrator Local Declarations
struct SubsurfaceOctreeNode {
    // SubsurfaceOctreeNode Methods
    SubsurfaceOctreeNode() {
        isLeaf = true;
        sumArea = 0.f;
        for (int i = 0; i < 8; ++i)
            ips[i] = NULL;
    }
    void Insert(const BBox &nodeBound, IrradiancePoint *ip,
                MemoryArena &arena) {
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        if (isLeaf) {
            // Add _IrradiancePoint_ to leaf octree node
            for (int i = 0; i < 8; ++i) {
                if (!ips[i]) {
                    ips[i] = ip;
                    return;
                }
            }

            // Convert leaf node to interior node, redistribute points
            isLeaf = false;
            IrradiancePoint *localIps[8];
            for (int i = 0; i < 8; ++i) {
                localIps[i] = ips[i];
                children[i] = NULL;
            }
            for (int i = 0; i < 8; ++i)  {
                IrradiancePoint *ip = localIps[i];
                // Add _IrradiancePoint_ _ip_ to interior octree node
                int child = (ip->p.x > pMid.x ? 4 : 0) +
                    (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
                if (!children[child])
                    children[child] = arena.Alloc<SubsurfaceOctreeNode>();
                BBox childBound = octreeChildBound(child, nodeBound, pMid);
                children[child]->Insert(childBound, ip, arena);
            }
            /* fall through to interior case to insert the new point... */
        }
        // Add _IrradiancePoint_ _ip_ to interior octree node
        int child = (ip->p.x > pMid.x ? 4 : 0) +
            (ip->p.y > pMid.y ? 2 : 0) + (ip->p.z > pMid.z ? 1 : 0);
        if (!children[child])
            children[child] = arena.Alloc<SubsurfaceOctreeNode>();
        BBox childBound = octreeChildBound(child, nodeBound, pMid);
        children[child]->Insert(childBound, ip, arena);
    }
    void InitHierarchy() {
        if (isLeaf) {
            // Init _SubsurfaceOctreeNode_ leaf from _IrradiancePoint_s
            float sumWt = 0.f;
            uint32_t i;
            for (i = 0; i < 8; ++i) {
                if (!ips[i]) break;
                float wt = ips[i]->E.y();
                E += ips[i]->E;
                p += wt * ips[i]->p;
                sumWt += wt;
                sumArea += ips[i]->area;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= i;
        }
        else {
            // Init interior _SubsurfaceOctreeNode_
            float sumWt = 0.f;
            uint32_t nChildren = 0;
            for (uint32_t i = 0; i < 8; ++i) {
                if (!children[i]) continue;
                ++nChildren;
                children[i]->InitHierarchy();
                float wt = children[i]->E.y();
                E += children[i]->E;
                p += wt * children[i]->p;
                sumWt += wt;
                sumArea += children[i]->sumArea;
            }
            if (sumWt > 0.f) p /= sumWt;
            E /= nChildren;
        }
    }
    Spectrum Mo(const BBox &nodeBound, const Point &pt, const float maxError, const vector<Spectrum> &R12, Spectrum &maxMo);

    // SubsurfaceOctreeNode Public Data
    Point p;
    bool isLeaf;
    Spectrum E;
    float sumArea;
    union {
        SubsurfaceOctreeNode *children[8];
        IrradiancePoint *ips[8];
    };
};


struct DiffusionReflectance {
    // DiffusionReflectance Public Methods
    DiffusionReflectance(const Spectrum &sigma_a, const Spectrum &sigmap_s,
                         float eta, float thickness) {
        A = (1.f + Fdr(eta)) / (1.f - Fdr(eta));
        sigmap_t = sigma_a + sigmap_s;
        sigma_tr = Sqrt(3.f * sigma_a * sigmap_t);
        alphap = sigmap_s / sigmap_t;
        zpos = Spectrum(1.f) / sigmap_t;
        zneg = -zpos * (1.f + (4.f/3.f) * A);
    }
    void operator()(float d2, vector<Spectrum> &Rd, vector<Spectrum>& Td) const {
        Spectrum dpos = Sqrt(Spectrum(d2) + zpos * zpos);
        Spectrum dneg = Sqrt(Spectrum(d2) + zneg * zneg);
        Spectrum R = (alphap / (4.f * M_PI)) *
            ((zpos * (dpos * sigma_tr + Spectrum(1.f)) *
              Exp(-sigma_tr * dpos)) / (dpos * dpos * dpos) -
             (zneg * (dneg * sigma_tr + Spectrum(1.f)) *
              Exp(-sigma_tr * dneg)) / (dneg * dneg * dneg));
        Rd.push_back(R.Clamp());
        Td.push_back(R.Clamp());
    }

    // DiffusionReflectance Data
    Spectrum zpos, zneg, sigmap_t, sigma_tr, alphap;
    float A;
};

struct MultipoleDR {
    // Initializes constants for Multipole DiffusionReflectance
    MultipoleDR(const Spectrum &sigma_a, const Spectrum &sigmap_s, float eta, float thickness) {
        A = (1.f + Fdr(eta)) / (1.f - Fdr(eta));
        sigmap_t = sigma_a + sigmap_s;
        sigma_tr = Sqrt(3.f * sigma_a * sigmap_t);        
        alphap = sigmap_s / sigmap_t;

        
        Spectrum lfree = Spectrum(1.f) / sigmap_t;
        depth = Spectrum(thickness);

        Spectrum zb = (2.f/3.f) * A * lfree;
        int shift = (NUM_DIPOLE/2);
        for (int i = 0; i < NUM_DIPOLE; i++) {
            int j = i-shift;
            zpos[i] = 2*j*(depth + 2*zb) + lfree;
            zneg[i] = 2*j*(depth + 2*zb) - lfree - 2*zb;
        }
    }

    // Computes R(r) for the given Multipole DiffusionReflectance constants
    void operator()(float d2, vector<Spectrum> &Rd, vector<Spectrum> &Td) const {
        Spectrum R = Spectrum(0.f);
        Spectrum T = Spectrum(0.f);
        for (int i = 0; i < NUM_DIPOLE; i++) {
            Spectrum dpos = Sqrt(Spectrum(d2) + zpos[i] * zpos[i]);
            Spectrum dneg = Sqrt(Spectrum(d2) + zneg[i] * zneg[i]);
            R += (alphap / (4.f * M_PI)) *
                ((zpos[i] * (dpos * sigma_tr + Spectrum(1.f)) *
                  Exp(-sigma_tr * dpos)) / (dpos * dpos * dpos) -
                 (zneg[i] * (dneg * sigma_tr + Spectrum(1.f)) *
                  Exp(-sigma_tr * dneg)) / (dneg * dneg * dneg));

            T += (alphap / (4.f * M_PI)) *
                (((depth - zpos[i]) * (dpos * sigma_tr + Spectrum(1.f)) *
                  Exp(-sigma_tr * dpos)) / (dpos * dpos * dpos) -
                 ((depth - zneg[i]) * (dneg * sigma_tr + Spectrum(1.f)) *
                  Exp(-sigma_tr * dneg)) / (dneg * dneg * dneg));            
        }

        Rd.push_back(R.Clamp());
        Td.push_back(T.Clamp());
    }

    // Multipole DiffusionReflectance Data
    Spectrum zpos[NUM_DIPOLE], zneg[NUM_DIPOLE], sigmap_t, sigma_tr, alphap, depth;
    float A;
};


// performs convolution
void convolve(vector<Spectrum> &s1, vector<Spectrum> &s2, vector<Spectrum> &result) {
    int size = s1.size() + s2.size() - 1;
    for (int i = 0; i < size; i++) {
        int i1 = i;
        result.push_back(Spectrum(0.f));
        for (int j = 0; j < int(s2.size()); j++) {
            if (i1 >= 0 && i1 < int(s1.size())) {
                result[i] += s1[i1] * s2[j];
            }
            i1--;
        }
        result[i] *= RADIUS_PREC;
    }
}

// DipoleSubsurfaceIntegrator Method Definitions
DipoleSubsurfaceIntegrator::~DipoleSubsurfaceIntegrator() {
    delete[] lightSampleOffsets;
    delete[] bsdfSampleOffsets;
}


void DipoleSubsurfaceIntegrator::RequestSamples(Sampler *sampler, Sample *sample,
        const Scene *scene) {
    // Allocate and request samples for sampling all lights
    uint32_t nLights = scene->lights.size();
    lightSampleOffsets = new LightSampleOffsets[nLights];
    bsdfSampleOffsets = new BSDFSampleOffsets[nLights];
    for (uint32_t i = 0; i < nLights; ++i) {
        const Light *light = scene->lights[i];
        int nSamples = light->nSamples;
        if (sampler) nSamples = sampler->RoundSize(nSamples);
        lightSampleOffsets[i] = LightSampleOffsets(nSamples, sample);
        bsdfSampleOffsets[i] = BSDFSampleOffsets(nSamples, sample);
    }
}


void DipoleSubsurfaceIntegrator::Preprocess(const Scene *scene,
        const Camera *camera, const Renderer *renderer) {
    if (scene->lights.size() == 0) return;
    vector<SurfacePoint> pts;
    // Get _SurfacePoint_s for translucent objects in scene
    if (filename != "") {
        // Initialize _SurfacePoint_s from file
        vector<float> fpts;
        if (ReadFloatFile(filename.c_str(), &fpts)) {
            if ((fpts.size() % 8) != 0)
                Error("Excess values (%d) in points file \"%s\"", int(fpts.size() % 8),
                      filename.c_str());
            for (u_int i = 0; i < fpts.size(); i += 8)
                pts.push_back(SurfacePoint(Point(fpts[i], fpts[i+1], fpts[i+2]),
                                           Normal(fpts[i+3], fpts[i+4], fpts[i+5]),
                                           fpts[i+6], fpts[i+7]));
        }
    }
    if (pts.size() == 0) {
        Point pCamera = camera->CameraToWorld(camera->shutterOpen,
                                              Point(0, 0, 0));
        FindPoissonPointDistribution(pCamera, camera->shutterOpen,
                                     minSampleDist, scene, &pts);
    }

    // Compute irradiance values at sample points
    camera_pr = camera;
    RNG rng;
    MemoryArena arena;
    PBRT_SUBSURFACE_STARTED_COMPUTING_IRRADIANCE_VALUES();
    ProgressReporter progress(pts.size(), "Computing Irradiances");
    for (uint32_t i = 0; i < pts.size(); ++i) {
        SurfacePoint &sp = pts[i];
        Spectrum E(0.f);
        for (uint32_t j = 0; j < scene->lights.size(); ++j) {
            // Add irradiance from light at point
            const Light *light = scene->lights[j];
            Spectrum Elight = 0.f;
            int nSamples = RoundUpPow2(light->nSamples);
            uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
            uint32_t compScramble = rng.RandomUInt();
            for (int s = 0; s < nSamples; ++s) {
                float lpos[2];
                Sample02(s, scramble, lpos);
                float lcomp = VanDerCorput(s, compScramble);
                LightSample ls(lpos[0], lpos[1], lcomp);
                Vector wi;
                float lightPdf;
                VisibilityTester visibility;
                Spectrum Li = light->Sample_L(sp.p, sp.rayEpsilon,
                    ls, camera->shutterOpen, &wi, &lightPdf, &visibility);
                if (Dot(wi, sp.n) <= 0.) continue;
                if (Li.IsBlack() || lightPdf == 0.f) continue;
                Li *= visibility.Transmittance(scene, renderer, NULL, rng, arena);
                if (visibility.Unoccluded(scene))
                    Elight += Li * AbsDot(wi, sp.n) / lightPdf;
            }
            E += Elight / nSamples;
        }
        irradiancePoints.push_back(IrradiancePoint(sp, E));
        PBRT_SUBSURFACE_COMPUTED_IRRADIANCE_AT_POINT(&sp, &E);
        arena.FreeAll();
        progress.Update();
    }
    progress.Done();
    PBRT_SUBSURFACE_FINISHED_COMPUTING_IRRADIANCE_VALUES();

    // Create octree of clustered irradiance samples
    octree = octreeArena.Alloc<SubsurfaceOctreeNode>();
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octreeBounds = Union(octreeBounds, irradiancePoints[i].p);
    for (uint32_t i = 0; i < irradiancePoints.size(); ++i)
        octree->Insert(octreeBounds, &irradiancePoints[i], octreeArena);
    octree->InitHierarchy(); 

    // Precompute reflectances and transmittances
    MultipoleDR Rd1(sigma_a_1, sigma_prime_s_1, eta_1, thickness_epi);
    MultipoleDR Rd2(sigma_a_2, sigma_prime_s_2, eta_2, thickness_derm);
    DiffusionReflectance Rd3(sigma_a_3, sigma_prime_s_3, eta_3, thickness_b_derm);

    vector<Spectrum> R1, R2, R3, T1, T2, T3;
    for (int i = 0; i < D12_SIZE; i++) {
        float dist = (float(i) * RADIUS_PREC) * (float(i) * RADIUS_PREC);
        Rd1(dist, R1, T1);
        Rd2(dist, R2, T2);
        Rd3(dist, R3, T3);
    }

    vector<Spectrum> R2T1, T1R2T1, R2R1, T1T2;
    convolve(R2, T1, R2T1);
    assert(R2T1.size() == R2.size() + T1.size() - 1);
    convolve(T1, R2T1, T1R2T1);
    assert(T1R2T1.size() == T1.size() + R2T1.size() - 1);
    convolve(R2, R1, R2R1);
    assert(R2R1.size() == R2.size() + R1.size() - 1);
    convolve(T1, T2, T1T2);

    vector<Spectrum> R12, T12;
    for (int i = 0; i < D12_SIZE; i++) {
        Spectrum temp = (T1R2T1[i]/(Spectrum(1.0f) - R2R1[i]));
        Spectrum temp2 = T1T2[i]/(Spectrum(1.0f) - R2R1[i]);
        R12.push_back(R1[i] + temp);
        T12.push_back(temp2);
        R1[i].Print();
        R2[i].Print();
        R2R1[i].Print();
        R12[i].Print();
        T12[i].Print();
        printf("\n");
    }

    vector<Spectrum> R3T12, T12R3T12, R3R12, T12T3;
    convolve(R3, T12, R3T12);
    assert(R3T12.size() == T12.size() + R3.size() -1);
    convolve(T12, R3T12, T12R3T12);
    assert(T12R3T12.size() == T12.size() + R3T12.size() -1);
    convolve(R3, R12, R3R12);
    assert(R3R12.size() == R3.size() + R12.size() - 1);
    convolve(T12, T3, T12T3);
    assert(T12T3.size() == T12.size() + T3.size() - 1);

    for (int i = 0; i < D12_SIZE; i++) {
        Spectrum temp = (T12R3T12[i]/(Spectrum(1.0f) - R3R12[i]));
        R3[i].Print();
        R12[i].Print();
        R3R12[i].Print();
        R123.push_back(R12[i] + temp);
        R123[i].Print();
        printf("\n");
    }    
    //exit(-1);

    // printf("R12\n");
    // for (int i = 0; i < D12_SIZE; i++) {
    //     R12[i].Print();
    // }
}


Spectrum DipoleSubsurfaceIntegrator::Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect,
        const Sample *sample, RNG &rng, MemoryArena &arena) const {  
    Spectrum L(0.);
    Vector wo = -ray.d;
    // Compute emitted light if ray hit an area light source
    L += isect.Le(wo);

    // Evaluate BSDF at hit point
    BSDF *bsdf = isect.GetBSDF(ray, arena);
    const Point &p = bsdf->dgShading.p;
    const Normal &n = bsdf->dgShading.nn;
    //Evaluate BSSRDF and possibly compute subsurface scattering
    BSSRDF *bssrdf = isect.GetBSSRDF(ray, arena);
    Spectrum temp = Spectrum(0.f);
    if (bssrdf && octree) {
        Spectrum sigma_a  = bssrdf->sigma_a();
        Spectrum sigmap_s = bssrdf->sigma_prime_s();
        Spectrum sigmap_t = sigmap_s + sigma_a;
        if (!sigmap_t.IsBlack()) {
            // Use hierarchical integration to evaluate reflection from dipole model
            PBRT_SUBSURFACE_STARTED_OCTREE_LOOKUP(const_cast<Point *>(&p));
            Spectrum maxMo = Spectrum(0.f);
            Spectrum Mo = octree->Mo(octreeBounds, p, maxError, R123, maxMo);
            Mo = Mo - maxMo;
            //Mo.Print();
            FresnelDielectric fresnel(1.f, bssrdf->eta());
            Spectrum Ft = Spectrum(1.f) - fresnel.Evaluate(AbsDot(wo, n));
            float Fdt = 1.f - Fdr(bssrdf->eta());

            temp = (INV_PI * Ft) * (Fdt * Mo);
            Spectrum E(0.f);
            for (uint32_t j = 0; j < scene->lights.size(); ++j) {
                // Add irradiance from light at point
                const Light *light = scene->lights[j];
                Spectrum Elight = 0.f;
                int nSamples = RoundUpPow2(light->nSamples);
                uint32_t scramble[2] = { rng.RandomUInt(), rng.RandomUInt() };
                uint32_t compScramble = rng.RandomUInt();
                for (int s = 0; s < nSamples; ++s) {
                    float lpos[2];
                    Sample02(s, scramble, lpos);
                    float lcomp = VanDerCorput(s, compScramble);
                    LightSample ls(lpos[0], lpos[1], lcomp);
                    Vector wi;
                    float lightPdf;
                    VisibilityTester visibility;
                    Spectrum Li = light->Sample_L(p, isect.rayEpsilon,
                        ls, camera_pr->shutterOpen, &wi, &lightPdf, &visibility);
                    if (Dot(wi, n) <= 0.) continue;
                    if (Li.IsBlack() || lightPdf == 0.f) continue;
                    Li *= visibility.Transmittance(scene, renderer, NULL, rng, arena);
                    if (visibility.Unoccluded(scene))
                        Elight += Li * AbsDot(wi, n) / lightPdf;
                }
                E += Elight / nSamples;
            }

            temp += E * R123[rand()%30+50] * INV_PI * Ft * Fdt;
            //R12[rand()%100 + 20].Print();
            //E.Print();
            //Mo.Print();
            //temp.Print();
            L += temp;
            PBRT_SUBSURFACE_FINISHED_OCTREE_LOOKUP();
        }
    }
    L += UniformSampleAllLights(scene, renderer, arena, p, n,
        wo, isect.rayEpsilon, ray.time, bsdf, sample, rng, lightSampleOffsets,
        bsdfSampleOffsets);
    if (ray.depth < maxSpecularDepth) {
        // Trace rays for specular reflection and refraction
        L += SpecularReflect(ray, bsdf, rng, isect, renderer, scene, sample,
                             arena);
        L += SpecularTransmit(ray, bsdf, rng, isect, renderer, scene, sample,
                              arena);
    }

    //Spectrum temp2 = (temp / L);
    //if (temp2.X() < 0.9 && temp2.Y() < 0.9 && temp2.Z() < 0.9 && temp2.X() > 0.1 && temp2.Y() > 0.1 && temp2.Z() > 0.1) return Spectrum(20.0f);
    // printf("temp : ");
    // temp.Print();
    // printf("L : ");
    // L.Print();

    //printf("temp / L : ");
    //temp2.Print();
    // if (temp2.X() > 1.f) {
    //     printf("temp : ");
    //     temp.Print();
    //     printf("L : ");
    //     L.Print();
    // }
    return L;
}


Spectrum SubsurfaceOctreeNode::Mo(const BBox &nodeBound, const Point &pt, const float maxError, const vector<Spectrum> &R12, Spectrum &maxMo) {
    // Compute $M_\roman{o}$ at node if error is low enough
    float dw = sumArea / DistanceSquared(pt, p);
    if (dw < maxError && !nodeBound.Inside(pt))
    {
        PBRT_SUBSURFACE_ADDED_INTERIOR_CONTRIBUTION(const_cast<SubsurfaceOctreeNode *>(this));
        int index = 0;
        float dist = sqrt(DistanceSquared(pt, p));
        if (dist >= MAX_RADIUS) {
            index = D12_SIZE - 1;
        /*} else if (dist <= RADIUS_PREC) {
            index = 0;*/
        } else {
            index = int(dist*100);
            if (index >= D12_SIZE) index = D12_SIZE-1;
        }
        Spectrum temp = R12[index] * E * sumArea;
        if (temp.X() > maxMo.X()) maxMo = temp;
        return temp;
    }

    // Otherwise compute $M_\roman{o}$ from points in leaf or recursively visit children
    Spectrum Mo = 0.f;
    if (isLeaf) {
        // Accumulate $M_\roman{o}$ from leaf node
        for (int i = 0; i < 8; ++i) {
            if (!ips[i]) break;
            PBRT_SUBSURFACE_ADDED_POINT_CONTRIBUTION(const_cast<IrradiancePoint *>(ips[i]));
            int index = 0;
            float dist = sqrt(DistanceSquared(pt, ips[i]->p));
            if (dist >= MAX_RADIUS) {
                index = D12_SIZE - 1;
            /*} else if (dist <= RADIUS_PREC) {
                index = 0;*/
            } else {
                index = int(dist*100);
                if (index >= D12_SIZE) index = D12_SIZE-1;
            }
            Spectrum temp = R12[index] * ips[i]->E * ips[i]->area;
            if (temp.X() > maxMo.X()) maxMo = temp;
            Mo += temp;
        }
    }
    else {
        // Recursively visit children nodes to compute $M_\roman{o}$
        Point pMid = .5f * nodeBound.pMin + .5f * nodeBound.pMax;
        for (int child = 0; child < 8; ++child) {
            if (!children[child]) continue;
            BBox childBound = octreeChildBound(child, nodeBound, pMid);
            Mo += children[child]->Mo(childBound, pt, maxError, R12, maxMo);
        }
    }
    return Mo;
}


DipoleSubsurfaceIntegrator *CreateDipoleSubsurfaceIntegrator(const ParamSet &params) {
    int maxDepth = params.FindOneInt("maxdepth", 5);
    float maxError = params.FindOneFloat("maxerror", .05f);
    float minDist = params.FindOneFloat("minsampledistance", .25f);
    float eta_1 = params.FindOneFloat("eta_1", 1.4f);
    float eta_2 = params.FindOneFloat("eta_2", 1.f);
    float eta_3 = params.FindOneFloat("eta_3", 1.f);
    float thickness_epi = params.FindOneFloat("thickness_epi", 0.25f);
    float thickness_derm = params.FindOneFloat("thickness_derm", 20.f);
    float thickness_b_derm = params.FindOneFloat("thickness_b_derm", 30.f);
    Spectrum sigma_a_1 = params.FindOneSpectrum("sigma_a_1", Spectrum(0.25f));
    Spectrum sigma_a_2 = params.FindOneSpectrum("sigma_a_2", Spectrum(0.25f));
    Spectrum sigma_a_3 = params.FindOneSpectrum("sigma_a_3", Spectrum(0.25f));
    Spectrum sigma_prime_s_1 = params.FindOneSpectrum("sigma_prime_s_1", Spectrum(0.25f));
    Spectrum sigma_prime_s_2 = params.FindOneSpectrum("sigma_prime_s_2", Spectrum(0.25f));
    Spectrum sigma_prime_s_3 = params.FindOneSpectrum("sigma_prime_s_3", Spectrum(0.25f));

    string pointsfile = params.FindOneFilename("pointsfile", "");
    if (PbrtOptions.quickRender) { maxError *= 4.f; minDist *= 4.f; }
    return new DipoleSubsurfaceIntegrator(maxDepth, maxError, minDist, 
        eta_1, eta_2, eta_3, thickness_epi, thickness_derm, thickness_b_derm, sigma_a_1, sigma_a_2, sigma_a_3,
        sigma_prime_s_1, sigma_prime_s_2, sigma_prime_s_3, pointsfile);
}


