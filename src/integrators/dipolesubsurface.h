
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

#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_INTEGRATORS_DIPOLESUBSURFACE_H
#define PBRT_INTEGRATORS_DIPOLESUBSURFACE_H

// integrators/dipolesubsurface.h*
#include "pbrt.h"
#include "integrator.h"
#include "kdtree.h"
#include "renderers/surfacepoints.h"
struct SubsurfaceOctreeNode;

// DipoleSubsurfaceIntegrator Helper Declarations
struct IrradiancePoint {
    IrradiancePoint() { }
    IrradiancePoint(const SurfacePoint &sp, const Spectrum &ee)
        : p(sp.p), n(sp.n), E(ee), area(sp.area),
          rayEpsilon(sp.rayEpsilon) { }
    Point p;
    Normal n;
    Spectrum E;
    float area, rayEpsilon;
};



// DipoleSubsurfaceIntegrator Declarations
class DipoleSubsurfaceIntegrator : public SurfaceIntegrator {
public:
    // DipoleSubsurfaceIntegrator Public Methods
    DipoleSubsurfaceIntegrator(int mdepth, float merror, float mindist, 
                                float e1, float e2, float e3, float t_epi, float t_derm, float t_b_derm, 
                                Spectrum s_a1, Spectrum s_a2, Spectrum s_a3, 
                                Spectrum s_p1, Spectrum s_p2, Spectrum s_p3, 
                                const string &fn) {
        maxSpecularDepth = mdepth;
        maxError = merror;
        minSampleDist = mindist;
        eta_1 = e1;
        eta_2 = e2;
        eta_3 = e3;
        thickness_epi = t_epi;
        thickness_derm = t_derm;
        thickness_b_derm = t_b_derm;
        sigma_a_1 = s_a1;
        sigma_a_2 = s_a2;
        sigma_a_3 = s_a3;
        sigma_prime_s_1 = s_p1;
        sigma_prime_s_2 = s_p2;
        sigma_prime_s_3 = s_p3;
        filename = fn;
        octree = NULL;
    }
    ~DipoleSubsurfaceIntegrator();
    Spectrum Li(const Scene *scene, const Renderer *renderer,
        const RayDifferential &ray, const Intersection &isect, const Sample *sample,
        RNG &rng, MemoryArena &arena) const;
    void RequestSamples(Sampler *sampler, Sample *sample, const Scene *scene);
    void Preprocess(const Scene *, const Camera *, const Renderer *);

private:
    // DipoleSubsurfaceIntegrator Private Data
    vector<Spectrum> R123;
    int maxSpecularDepth;
    float maxError, minSampleDist, eta_1, eta_2, eta_3, thickness_epi, thickness_derm, thickness_b_derm;
    Spectrum sigma_a_1, sigma_a_2, sigma_a_3, sigma_prime_s_1, sigma_prime_s_2, sigma_prime_s_3;
    string filename;
    vector<IrradiancePoint> irradiancePoints;
    BBox octreeBounds;
    SubsurfaceOctreeNode *octree;
    MemoryArena octreeArena;
    const Camera *camera_pr;

    // Declare sample parameters for light source sampling
    LightSampleOffsets *lightSampleOffsets;
    BSDFSampleOffsets *bsdfSampleOffsets;
};


DipoleSubsurfaceIntegrator *CreateDipoleSubsurfaceIntegrator(const ParamSet &params);

#endif // PBRT_INTEGRATORS_DIPOLESUBSURFACE_H
