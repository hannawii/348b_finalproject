// materials/skin.cpp*
#include "stdafx.h"
#include "materials/skin.h"
#include "spectrum.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

// UberMaterial Method Definitions
BSDF *SkinMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    // Allocate _BSDF_, possibly doing bump mapping with _bumpMap_
    // DifferentialGeometry dgs;
    // if (bumpMap)
    //     Bump(bumpMap, dgGeom, dgShading, &dgs);
    // else
    //     dgs = dgShading;
    // BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);

    // Spectrum op = opacity->Evaluate(dgs).Clamp();
    // if (op != Spectrum(1.)) {
    //     BxDF *tr = BSDF_ALLOC(arena, SpecularTransmission)(-op + Spectrum(1.), 1., 1.);
    //     bsdf->Add(tr);
    // }

    // Spectrum kd = op * Kd->Evaluate(dgs).Clamp();
    // if (!kd.IsBlack()) {
    //     BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
    //     bsdf->Add(diff);
    // }

    // float e = eta->Evaluate(dgs);
    // Spectrum ks = op * Ks->Evaluate(dgs).Clamp();
    // if (!ks.IsBlack()) {
    //     Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(e, 1.f);
    //     float rough = roughness->Evaluate(dgs);
    //     BxDF *spec = BSDF_ALLOC(arena, Microfacet)(ks, fresnel, BSDF_ALLOC(arena, Blinn)(1.f / rough));
    //     bsdf->Add(spec);
    // }

    // Spectrum kr = op * Kr->Evaluate(dgs).Clamp();
    // if (!kr.IsBlack()) {
    //     Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(e, 1.f);
    //     bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(kr, fresnel));
    // }

    // Spectrum kt = op * Kt->Evaluate(dgs).Clamp();
    // if (!kt.IsBlack())
    //     bsdf->Add(BSDF_ALLOC(arena, SpecularTransmission)(kt, e, 1.f));

    // return bsdf;
}

BSSRDF *SkinMaterial::GetBSSRDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    // float e = eta->Evaluate(dgShading);
    // float mfp = meanfreepath->Evaluate(dgShading);
    // Spectrum kd = Kd->Evaluate(dgShading).Clamp();
    // Spectrum sigma_a, sigma_prime_s;
    // SubsurfaceFromDiffuse(kd, mfp, e, &sigma_a, &sigma_prime_s);
    // return BSDF_ALLOC(arena, BSSRDF)(sigma_a, sigma_prime_s, e);
}


SkinMaterial *CreateSkinMaterial(const Transform &xform, const TextureParams &mp) {
    // Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    // Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    // Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(0.f));
    // Reference<Texture<Spectrum> > Kt = mp.GetSpectrumTexture("Kt", Spectrum(0.f));
    // Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", .1f);
    // Reference<Texture<float> > eta = mp.GetFloatTexture("index", 1.5f);
    // Reference<Texture<Spectrum> > opacity = mp.GetSpectrumTexture("opacity", 1.f);
    // Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");
    // return new SkinMaterial(Kd, Ks, Kr, Kt, roughness, opacity, eta, bumpMap);
}