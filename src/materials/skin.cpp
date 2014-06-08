// materials/skin.cpp*
#include "stdafx.h"
#include "materials/skin.h"
#include "spectrum.h"
#include "volume.h"
#include "reflection.h"
#include "texture.h"
#include "paramset.h"

// SkinMaterial Method Definitions
BSDF *SkinMaterial::GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    //Allocate _BSDF_, possibly doing bump mapping with _bumpMap_

    // Spectrum R = Kr->Evaluate(dgs).Clamp();
    // float e = eta->Evaluate(dgs);
    // if (!R.IsBlack())
    //     bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(R,
    //         BSDF_ALLOC(arena, FresnelDielectric)(1., e)));

    DifferentialGeometry dgs;
    if (bumpMap)
        Bump(bumpMap, dgGeom, dgShading, &dgs);
    else
        dgs = dgShading;
    BSDF *bsdf = BSDF_ALLOC(arena, BSDF)(dgs, dgGeom.nn);

    Spectrum op = opacity->Evaluate(dgs).Clamp();
    // if (op != Spectrum(1.)) {
    //     BxDF *tr = BSDF_ALLOC(arena, SpecularTransmission)(-op + Spectrum(1.), 1., 1.);
    //     bsdf->Add(tr);
    // }

    Spectrum kd = Kd->Evaluate(dgs).Clamp();
    if (!kd.IsBlack()) {
        BxDF *diff = BSDF_ALLOC(arena, Lambertian)(kd);
        bsdf->Add(diff);
    }

    float e = eta->Evaluate(dgs);
    Spectrum ks = Ks->Evaluate(dgs).Clamp();
    if (!ks.IsBlack()) {
        Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(e, 1.f);
        float rough = roughness->Evaluate(dgs);
        // bsdf->Add(BSDF_ALLOC(arena, FresnelBlend)(kd, ks, BSDF_ALLOC(arena, Blinn)(4.f)));
        BxDF *spec = BSDF_ALLOC(arena, Microfacet)(ks, fresnel, BSDF_ALLOC(arena, Blinn)(1.f / rough));
        bsdf->Add(spec);
    }

    // Spectrum kr = op * Kr->Evaluate(dgs).Clamp();
    // if (!kr.IsBlack()) {
    //     Fresnel *fresnel = BSDF_ALLOC(arena, FresnelDielectric)(e, 1.f);
    //     bsdf->Add(BSDF_ALLOC(arena, SpecularReflection)(kr, fresnel));
    //}

    // Spectrum kt = op * Kt->Evaluate(dgs).Clamp();
    // if (!kt.IsBlack())
    //     bsdf->Add(BSDF_ALLOC(arena, SpecularTransmission)(kt, e, 1.f));

    return bsdf;
}

BSSRDF *SkinMaterial::GetBSSRDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const {
    float e = eta->Evaluate(dgShading);
    return BSDF_ALLOC(arena, BSSRDF)(sigma_a->Evaluate(dgShading),
        sigma_prime_s->Evaluate(dgShading), e);
}


SkinMaterial *CreateSkinMaterial(const Transform &xform, const TextureParams &mp) {
    float sa_rgb[3] = { .0011f, .0024f, .014f }, sps_rgb[3] = { 2.55f, 3.21f, 3.77f };
    Spectrum sa = Spectrum::FromRGB(sa_rgb), sps = Spectrum::FromRGB(sps_rgb);
    string name = mp.FindString("name");
    bool found = GetVolumeScatteringProperties(name, &sa, &sps);
    if (name != "" && !found)
        Warning("Named material \"%s\" not found.  Using defaults.", name.c_str());
    float scale = mp.FindFloat("scale", 1.f);

    Reference<Texture<Spectrum> > sigma_a, sigma_prime_s;
    sigma_a = mp.GetSpectrumTexture("sigma_a", sa);
    sigma_prime_s = mp.GetSpectrumTexture("sigma_prime_s", sps);

    Reference<Texture<Spectrum> > Kd = mp.GetSpectrumTexture("Kd", Spectrum(0.25f));
    Reference<Texture<Spectrum> > Ks = mp.GetSpectrumTexture("Ks", Spectrum(0.25f));
    Reference<Texture<Spectrum> > Kr = mp.GetSpectrumTexture("Kr", Spectrum(0.f));
    Reference<Texture<Spectrum> > Kt = mp.GetSpectrumTexture("Kt", Spectrum(0.f));
    Reference<Texture<float> > roughness = mp.GetFloatTexture("roughness", 3.5f);
    Reference<Texture<float> > eta = mp.GetFloatTexture("index", 1.5f);
    Reference<Texture<Spectrum> > opacity = mp.GetSpectrumTexture("opacity", 1.f);
    Reference<Texture<float> > bumpMap = mp.GetFloatTextureOrNull("bumpmap");

    return new SkinMaterial(Kd, Ks, Kr, Kt, roughness, opacity, eta, bumpMap, sigma_a, sigma_prime_s, scale);
}