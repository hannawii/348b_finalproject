#if defined(_MSC_VER)
#pragma once
#endif

#ifndef PBRT_MATERIALS_SKIN_H
#define PBRT_MATERIALS_SKIN_H

// materials/skin.h*
#include "pbrt.h"
#include "material.h"

// SkinMaterial Declarations
class SkinMaterial : public Material {
public:
    SkinMaterial(Reference<Texture<Spectrum> > kd,
        Reference<Texture<Spectrum> > ks,
        Reference<Texture<Spectrum> > kr,
        Reference<Texture<Spectrum> > kt,
        Reference<Texture<float> > rough,
        Reference<Texture<Spectrum> > op,
        Reference<Texture<float> > e,
        Reference<Texture<float> > bump,
        Reference<Texture<Spectrum> > sa,
        Reference<Texture<Spectrum> > sps,
        float sc) {

        Kd = kd;
        Ks = ks;
        Kr = kr;
        Kt = kt;
        roughness = rough;
        opacity = op;
        eta = e;
        bumpMap = bump;
        scale = sc;
        sigma_a = sa;
        sigma_prime_s = sps;
    }
    BSDF *GetBSDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const;

    BSSRDF *GetBSSRDF(const DifferentialGeometry &dgGeom, const DifferentialGeometry &dgShading, MemoryArena &arena) const;
private:
    // SkinMaterial Private Data
    Reference<Texture<Spectrum> > Kd, Ks, Kr, Kt, opacity, sigma_a, sigma_prime_s;
    Reference<Texture<float> > roughness, eta, bumpMap;
    float scale;
};


SkinMaterial *CreateSkinMaterial(const Transform &xform,
        const TextureParams &mp);

#endif // PBRT_MATERIALS_SKIN_H