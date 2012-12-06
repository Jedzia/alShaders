#include <ai.h>
#include <cstring>
#include <iostream>
#include <ai_sampler.h>
#include <OpenEXR/ImathVec.h>
#include <OpenEXR/ImathMatrix.h>
#include <OpenEXR/ImathMatrixAlgo.h>

#include "alUtil.h"
#include "MIS.h"
#include "BeckmannMicrofacet.h"
#include "alSurface.h"
#include "Shadows.h"

AI_SHADER_NODE_EXPORT_METHODS(alSurfaceMtd)

#define GlossyMISBRDF AiCookTorranceMISBRDF
#define GlossyMISPDF AiCookTorranceMISPDF
#define GlossyMISSample AiCookTorranceMISSample

#define GlossyMISBRDF_wrap AiCookTorranceMISBRDF_wrap
#define GlossyMISPDF_wrap AiCookTorranceMISPDF_wrap
#define GlossyMISSample_wrap AiCookTorranceMISSample_wrap

#define GlossyMISCreateData AiCookTorranceMISCreateData

enum alSurfaceParams
{
	// diffuse
	p_diffuseScale=0,
	p_diffuseColor,
	p_diffuseRoughness,
	p_emissionScale,
	p_emissionColor,

	// sss
	p_sssMix,
	p_sssRadius,
	p_sssRadiusColor,
	p_sssScale,

	p_ssScale,
	p_ssTargetDistance,
	p_ssTargetColor,
	p_ssSpecifyCoefficients,
	p_ssScattering,
	p_ssAbsorption,
	p_ssUnitScale,
	p_ssDirection,

	p_diffuseExtraSamples,
	p_diffuseEnableCaustics,

	// specular
	p_specular1Scale,
	p_specular1Color,
	p_specular1Roughness,
	p_specular1Ior,
	p_specular1RoughnessDepthScale,
	p_specular1ExtraSamples,
	p_specular2Scale,
	p_specular2Color,
	p_specular2Roughness,
	p_specular2Ior,
	p_specular2RoughnessDepthScale,
	p_specular2ExtraSamples,

	// transmission
	p_transmissionScale,
	p_transmissionColor,
	p_transmissionLinkToSpecular1,
	p_transmissionRoughness,
	p_transmissionIor,
	p_transmissionRoughnessDepthScale,
	p_transmissionEnableCaustics,
	p_transmissionExtraSamples,
	p_absorptionEnable,
	p_absorptionDensity,
	p_absorptionColor,

	p_bump
};

node_parameters
{
	AiParameterFLT("diffuseScale", 1.0f );
	AiParameterRGB("diffuseColor", 0.18f, 0.18f, 0.18f );
	AiParameterFLT("diffuseRoughness", 0.0f );
	AiParameterFLT("emissionScale", 0.0f );
	AiParameterRGB("emissionColor", 1.0f, 1.0f, 1.0f);

	AiParameterFLT("sssMix", 0.0f );
	AiParameterFLT("sssRadius", 3.6f );
	AiParameterRGB("sssRadiusColor", .439f, .156f, .078f );
	AiMetaDataSetBool(mds, "sssRadiusColor", "always_linear", true);  // no inverse-gamma correction
	AiParameterFLT("sssScale", 1.0f );

	AiParameterFLT("ssScale", 0.0f );
	AiParameterFLT("ssTargetDistance", 3.6f);
	AiParameterRGB("ssTargetColor", .439f, .156f, .078f);
	AiParameterBOOL("ssSpecifyCoefficients", false);
	AiParameterRGB("ssScattering", 1.0f, 1.0f, 1.0f);
	AiParameterRGB("ssAbsorption", 1.0f, 1.0f, 1.0f);
	AiParameterFLT("ssUnitScale", 1.0f);
	AiParameterFLT("ssDirection", 0.0f);

	AiParameterINT("diffuseExtraSamples", 0);
	AiParameterBOOL("diffuseEnableCaustics", false);

	AiParameterFLT("specular1Scale", 1.0f );
	AiParameterRGB("specular1Color", 1.0f, 1.0f, 1.0f );
	AiParameterFLT("specular1Roughness", 0.3f );
	AiParameterFLT("specular1Ior", 1.4f );
	AiParameterFLT("specular1RoughnessDepthScale", 1.0f);
	AiParameterINT("specular1ExtraSamples", 0);

	AiParameterFLT("specular2Scale", 0.0f );
	AiParameterRGB("specular2Color", 1.0f, 1.0f, 1.0f );
	AiParameterFLT("specular2Roughness", 0.3f );
	AiParameterFLT("specular2Ior", 1.4f );
	AiParameterFLT("specular2RoughnessDepthScale", 1.0f);
	AiParameterINT("specular2ExtraSamples", 0);

	AiParameterFLT("transmissionScale", 0.0f );
	AiParameterRGB("transmissionColor", 1.0f, 1.0f, 1.0f );
	AiParameterBOOL("transmissionLinkToSpecular1", true);
	AiParameterFLT("transmissionRoughness", 0.1f );
	AiParameterFLT("transmissionIor", 1.4f );
	AiParameterFLT("transmissionRoughnessDepthScale", 1.0f);
	AiParameterBOOL("transmissionEnableCaustics", true);
	AiParameterINT("transmissionExtraSamples", 0);
	AiParameterBOOL("absorptionEnable", false);
	AiParameterFLT("absorptionDensity", 1.0f);
	AiParameterRGB("absorptionColor", 1.0f, 1.0f, 1.0f);
}


node_loader
{
	if (i>0) return 0;
	node->methods     = alSurfaceMtd;
	node->output_type = AI_TYPE_RGB;
	node->name        = "alSurface";
	node->node_type   = AI_NODE_SHADER;
	strcpy(node->version, AI_VERSION);
	return TRUE;
}

node_initialize
{
	ShaderData *data = (ShaderData*) AiMalloc(sizeof(ShaderData));
	AiNodeSetLocalData(node,data);
	data->diffuse_sampler = NULL;
	data->glossy_sampler = NULL;
	data->glossy2_sampler = NULL;
	data->refraction_sampler = NULL;
};

node_finish
{
	if (AiNodeGetLocalData(node))
	{
		ShaderData* data = (ShaderData*) AiNodeGetLocalData(node);

		AiSamplerDestroy(data->diffuse_sampler);
		AiSamplerDestroy(data->glossy_sampler);
		AiSamplerDestroy(data->glossy2_sampler);
		AiSamplerDestroy(data->refraction_sampler);

		AiFree((void*) data);
		AiNodeSetLocalData(node, NULL);
		}
}


node_update
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);
	AtNode *options   = AiUniverseGetOptions();
	data->GI_diffuse_depth = AiNodeGetInt(options, "GI_diffuse_depth");
	data->GI_reflection_depth = AiNodeGetInt(options, "GI_reflection_depth");
	data->GI_refraction_depth = AiNodeGetInt(options, "GI_refraction_depth");
	data->GI_glossy_depth = AiNodeGetInt(options, "GI_glossy_depth");
	data->GI_glossy_samples = AiNodeGetInt(options, "GI_glossy_samples");
	data->GI_diffuse_samples = AiNodeGetInt(options, "GI_diffuse_samples");

	// setup samples
	AiSamplerDestroy(data->diffuse_sampler);
	AiSamplerDestroy(data->glossy_sampler);
	AiSamplerDestroy(data->glossy2_sampler);
	AiSamplerDestroy(data->refraction_sampler);
	data->diffuse_sampler = AiSampler(data->GI_diffuse_samples+params[p_diffuseExtraSamples].INT, 2);
	data->glossy_sampler = AiSampler(data->GI_glossy_samples+params[p_specular1ExtraSamples].INT, 2);
	data->glossy2_sampler = AiSampler(data->GI_glossy_samples+params[p_specular2ExtraSamples].INT, 2);
	data->refraction_sampler = AiSampler(AiNodeGetInt(options, "GI_refraction_samples")+params[p_transmissionExtraSamples].INT, 2);
};


shader_evaluate
{
	ShaderData *data = (ShaderData*)AiNodeGetLocalData(node);

	// Evaluate bump;
	AtVector N_orig;
	AtVector Nf_orig;
	AtRGB bump = AiShaderEvalParamRGB( p_bump );

	// Initialize parameter temporaries
	AtRGB diffuseColor = AiShaderEvalParamRGB( p_diffuseColor ) * AiShaderEvalParamFlt( p_diffuseScale );
	AtFloat diffuseRoughness = AiShaderEvalParamFlt(p_diffuseRoughness);
	bool diffuseEnableCaustics = AiShaderEvalParamFlt(p_diffuseEnableCaustics);
	AtRGB emissionColor = AiShaderEvalParamRGB(p_emissionColor) * AiShaderEvalParamFlt(p_emissionScale);
	AtFloat sssMix = AiShaderEvalParamFlt( p_sssMix );
	AtRGB sssRadiusColor = AiShaderEvalParamRGB( p_sssRadiusColor );
	AtFloat sssRadius = AiShaderEvalParamFlt( p_sssRadius );
	AtFloat sssScale = AiShaderEvalParamFlt( p_sssScale );
	AtRGB specular1Color = AiShaderEvalParamRGB( p_specular1Color ) * AiShaderEvalParamFlt( p_specular1Scale );
	AtRGB specular2Color = AiShaderEvalParamRGB( p_specular2Color ) * AiShaderEvalParamFlt( p_specular2Scale );
	AtFloat roughness = AiShaderEvalParamFlt( p_specular1Roughness );
	roughness *= roughness;
	AtFloat roughness2 = AiShaderEvalParamFlt( p_specular2Roughness );
	roughness2 *= roughness2;
	AtFloat ior = AiShaderEvalParamFlt( p_specular1Ior );
	AtFloat eta = 1.0f / ior;
	AtFloat ior2 = AiShaderEvalParamFlt( p_specular2Ior );
	AtFloat eta2 = 1.0f / ior2;

	AtRGB ssScattering = AiShaderEvalParamRGB(p_ssScattering);
	AtRGB ssAbsorption = AiShaderEvalParamRGB(p_ssAbsorption);
	AtFloat ssUnitScale = AiShaderEvalParamFlt( p_ssUnitScale );
	AtFloat ssScale = AiShaderEvalParamFlt( p_ssScale );
	AtFloat ssDirection = AiShaderEvalParamFlt(p_ssDirection);
	AtFloat ssTargetDistance = AiShaderEvalParamFlt(p_ssTargetDistance);
	AtRGB ssTargetColor = AiShaderEvalParamRGB(p_ssTargetColor);
	bool ssSpecifyCoefficients = AiShaderEvalParamBool(p_ssSpecifyCoefficients);

	// precalculate scattering coefficients as we'll need them for shadows etc.
	AtRGB sigma_t = AI_RGB_BLACK;
	AtRGB sigma_s = AI_RGB_BLACK;
	AtRGB sigma_a = AI_RGB_BLACK;
	if (ssScale > IMPORTANCE_EPS)
	{
		if (ssSpecifyCoefficients)
		{
			sigma_s = ssScattering * ssUnitScale;
			sigma_a = ssAbsorption * ssUnitScale;
			sigma_t = sigma_s + sigma_a;
		}
		else
		{
			// do alpha inversion to construct scattering parameters
			AtRGB sigma_s_prime;
			// magic number of 0.125f comes from trying to estimate a matching mean free path to Arnold's cubic
			// scattering kernel, since 0.5 = (x/d)^-3
			alphaInversion(ssTargetColor*ssTargetDistance*0.125f, ssTargetDistance*0.125f, sigma_s_prime, sigma_a);
			sigma_s = sigma_s_prime / (1.0f - ssDirection);
		}
	}

	AtRGB transmissionColor = AiShaderEvalParamRGB(p_transmissionColor) * AiShaderEvalParamFlt(p_transmissionScale);

	AtFloat specular1RoughnessDepthScale = AiShaderEvalParamFlt(p_specular1RoughnessDepthScale);
	AtFloat specular2RoughnessDepthScale = AiShaderEvalParamFlt(p_specular2RoughnessDepthScale);
	AtFloat transmissionRoughnessDepthScale = AiShaderEvalParamFlt(p_transmissionRoughnessDepthScale);

	AtFloat transmissionRoughness;
	AtFloat transmissionIor;
	bool transmissionLinkToSpecular1 = AiShaderEvalParamBool(p_transmissionLinkToSpecular1);
	if (transmissionLinkToSpecular1)
	{
		transmissionRoughness = roughness;
		transmissionIor = ior;
	}
	else
	{
		transmissionRoughness = AiShaderEvalParamFlt(p_transmissionRoughness);
		transmissionRoughness *= transmissionRoughness;
		transmissionIor = AiShaderEvalParamFlt(p_transmissionIor);
	}
	bool transmissionEnableCaustics = AiShaderEvalParamBool(p_transmissionEnableCaustics);

	AtRGB absorption = AI_RGB_BLACK;
	if (AiShaderEvalParamBool(p_absorptionEnable))
	{
		absorption = (AI_RGB_WHITE - AiShaderEvalParamRGB(p_absorptionColor));
		absorption = max(rgb(AI_EPSILON), absorption) * AiShaderEvalParamFlt(p_absorptionDensity);
	}

	if (sg->Rt & AI_RAY_SHADOW)
	{
		// if the object is transmissive and
		AtRGB outOpacity = AI_RGB_WHITE;
		if (maxh(transmissionColor))
		{
			// check transmission through the surface
			AtFloat costheta = AiV3Dot(sg->Nf, -sg->Rd);
			AtFloat kt = 1.0f - fresnel(costheta, 1.0f/transmissionIor);
			if (kt >= IMPORTANCE_EPS) // else surface is fully reflective
			{
				if (maxh(sigma_t) > 0.0f)
				{
					AtPoint alsPreviousIntersection;
					AtRGB als_sigma_t = sigma_t;
					if (AiStateGetMsgPnt("alsPreviousIntersection", &alsPreviousIntersection))
					{
						AiStateGetMsgRGB("alsPrevious_sigma_t", &als_sigma_t);
						bool doExtinction = false;
						if (AiV3Dot(sg->N, sg->Rd) < 0.0f)
						{
							// ray is entering a closed volume
							bool alsInside;
							AiStateGetMsgBool("alsInside", &alsInside);
							if (alsInside)
							{
								// ray is entering an embedded volume
								doExtinction = true;
							}
							else
							{
								// shouldn't get here
							}

						}
						else
						{
							// ray is exiting a closed volume
							doExtinction = true;
						}

						if (doExtinction)
						{
							AtFloat z = AiV3Dist(sg->P, alsPreviousIntersection);
							outOpacity.r = expf(-z * als_sigma_t.r);
							outOpacity.g = expf(-z * als_sigma_t.g);
							outOpacity.b = expf(-z * als_sigma_t.b);
							outOpacity = 1.0f - (outOpacity*kt);
						}

					}
					else
					{
						// first intersection
						// tell the next shader invocation that we're now inside the surface and what our extinction
						// coefficient is
						AiStateSetMsgRGB("alsPrevious_sigma_t", sigma_t);
						AiStateSetMsgBool("alsInside", true);
					}
				}
				else // no extinction, shadows are fresnel only.
				{
					AiStateSetMsgRGB("alsPrevious_sigma_t", AI_RGB_BLACK);
					outOpacity = 1.0f - kt;
				}
			}
		}

		// store intersection position
		AiStateSetMsgPnt("alsPreviousIntersection", sg->P);
		sg->out_opacity = outOpacity;
		return;
	}

	// clamp roughnesses
	// TODO: fall back to single-ray solution when roughness is 0
	roughness = std::max(0.0001f, roughness);
	roughness2 = std::max(0.0001f, roughness2);
	transmissionRoughness = std::max(0.0001f, transmissionRoughness);

	// Initialize result temporaries
	AtRGB result_diffuseDirect = AI_RGB_BLACK;
	AtRGB result_glossyDirect = AI_RGB_BLACK;
	AtRGB result_glossy2Direct = AI_RGB_BLACK;
	AtRGB result_diffuseIndirect = AI_RGB_BLACK;
	AtRGB result_glossyIndirect = AI_RGB_BLACK;
	AtRGB result_glossy2Indirect = AI_RGB_BLACK;
	AtRGB result_sss = AI_RGB_BLACK;
	AtRGB result_ss = AI_RGB_BLACK;
	AtColor	result_transmission = AI_RGB_BLACK;
	AtColor result_emission = AI_RGB_BLACK;
	// Set up flags to early out of calculations based on where we are in the ray tree
	bool do_diffuse = true;
	bool do_glossy = true;
	bool do_glossy2 = true;
	bool do_ss = true;
	bool do_sss = true;
	bool do_transmission = true;
	AtInt glossy_samples = data->GI_glossy_samples;
	AtInt diffuse_samples = data->GI_diffuse_samples;

	if ( sg->Rr_diff > data->GI_diffuse_depth || maxh(diffuseColor) < IMPORTANCE_EPS)
	{
		do_diffuse = false;
	}


	if (sg->Rr_gloss > data->GI_glossy_depth
				|| (sg->Rr_diff > 0)			// disable glossy->diffuse caustics
				|| maxh(specular1Color) < IMPORTANCE_EPS				// skip evaluations that aren't important
				|| (sg->Rr_refr > 1 && !transmissionEnableCaustics))	// disable glossy->transmitted caustics
	{
		do_glossy = false;
	}

	if (sg->Rr_gloss > data->GI_glossy_depth
			|| (sg->Rr_diff > 0)			// disable glossy->diffuse caustics
			|| maxh(specular2Color) < IMPORTANCE_EPS				// skip evaluations that aren't important
			|| (sg->Rr_refr > 1 && !transmissionEnableCaustics))	// disable glossy->transmitted caustics
	{
		do_glossy2 = false;
	}

	if (sg->Rr_gloss > 0 && do_glossy)
	{
		roughness *= powf(specular1RoughnessDepthScale, sg->Rr_gloss);
		roughness2 *= powf(specular2RoughnessDepthScale, sg->Rr_gloss);
	}
	if (sg->Rr_refr > 0 && do_transmission)
	{
		transmissionRoughness *= powf(transmissionRoughnessDepthScale, sg->Rr_refr);
	}

	if ( sg->Rr_diff > 0 || sg->Rr_gloss > 1 || sssMix < 0.01f || !do_diffuse )
	{
		do_sss = false;
		sssMix = 0.0f;
	}

	if ((sg->Rr_diff > 0) || maxh(transmissionColor) < IMPORTANCE_EPS)
	{
		do_transmission = false;
	}

	// build a local frame for sampling
	AtVector U, V;
	AiBuildLocalFramePolar(&U, &V, &sg->N);

	AtVector wo = -sg->Rd;

	// Begin illumination calculation
	if (do_diffuse || do_glossy || do_glossy2)
	{
		// Create the BRDF data structures for MIS
		void* mis;
		mis = GlossyMISCreateData(sg,&U,&V,roughness,roughness);
		BrdfData_wrap brdfw;
		brdfw.brdf_data = mis;
		brdfw.sg = sg;
		brdfw.eta = eta;
		brdfw.V = wo;
		brdfw.N = sg->N;

		void* mis2;
		mis2 = GlossyMISCreateData(sg,&U,&V,roughness2,roughness2);
		BrdfData_wrap brdfw2;
		brdfw2.brdf_data = mis2;
		brdfw2.sg = sg;
		brdfw2.eta = eta2;
		brdfw2.V = wo;
		brdfw2.N = sg->N;

		void* dmis;
		dmis = AiOrenNayarMISCreateData(sg, diffuseRoughness);
		BrdfData_wrap brdfd;
		brdfd.brdf_data = dmis;
		brdfd.sg = sg;
		brdfd.eta = eta;
		brdfd.V = wo;
		brdfd.N = sg->N;

		// Light loop
		AiLightsPrepare(sg);
		while(AiLightsGetSample(sg))
		{
			if (do_diffuse)
			{
				result_diffuseDirect +=
				AiEvaluateLightSample(sg,&brdfd,AiOrenNayarMISSample_wrap,AiOrenNayarMISBRDF_wrap, AiOrenNayarMISPDF_wrap);
			}
			if (do_glossy)
			{
				result_glossyDirect +=
				AiEvaluateLightSample(sg,&brdfw,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap);
			}
			if (do_glossy2)
			{
				result_glossy2Direct +=
				AiEvaluateLightSample(sg,&brdfw2,GlossyMISSample_wrap,GlossyMISBRDF_wrap,GlossyMISPDF_wrap);
			}
		}

		// Multiply by the colors
		result_diffuseDirect *= diffuseColor;
		result_glossyDirect *= specular1Color;
		result_glossy2Direct *= specular2Color;

		// Sample BRDFS
		double samples[2];
		AtRay wi_ray;
		AtVector wi;
		AtScrSample scrs;
		AtVector H;
		float kr=1, kt=1;

		glossy_samples *= glossy_samples;
		if (do_glossy)
		{
			AtSamplerIterator* sampit = AiSamplerIterator(data->glossy_sampler, sg);
			AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
			AtInt count=0;
			while(AiSamplerGetSample(sampit, samples))
			{
				wi = GlossyMISSample(mis, samples[0], samples[1]);
				if (AiV3Dot(wi,sg->Nf) > 0.0f)
				{
					// get half-angle vector for fresnel
					wi_ray.dir = wi;
					AiV3Normalize(H, wi+brdfw.V);
					kr = fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta);
					if (kr > IMPORTANCE_EPS) // only trace a ray if it's going to matter
					{
						AiTrace(&wi_ray, &scrs);
						result_glossyIndirect +=
						scrs.color*GlossyMISBRDF(mis, &wi) / GlossyMISPDF(mis, &wi) * kr;
					}
				}
				count++;
			}
			if (count) result_glossyIndirect /= float(count);
			result_glossyIndirect *= specular1Color;
		} // if (do_glossy)

		if (do_glossy2)
		{
			AtSamplerIterator* sampit = AiSamplerIterator(data->glossy2_sampler, sg);
			AiMakeRay(&wi_ray, AI_RAY_GLOSSY, &sg->P, NULL, AI_BIG, sg);
			AtInt count=0;
			while(AiSamplerGetSample(sampit, samples))
			{
				wi = GlossyMISSample(mis2, samples[0], samples[1]);
				if (AiV3Dot(wi,sg->Nf) > 0.0f)
				{
					wi_ray.dir = wi;
					AiV3Normalize(H, wi+brdfw2.V);
					// add the fresnel for this layer
					kr = fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta2);
					// attenuate by the fresnel from the layer above
					kr *= 1.0f - fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta);
					if (kr > IMPORTANCE_EPS) // only trace a ray if it's going to matter
					{
						AiTrace(&wi_ray, &scrs);
						result_glossy2Indirect +=
						scrs.color*GlossyMISBRDF(mis2, &wi) / GlossyMISPDF(mis2, &wi) * kr;
					}
				}
				count++;
			}
			if (count) result_glossy2Indirect /= float(count);
			result_glossy2Indirect *= specular2Color;
		} // if (do_glossy2)

		if ( do_diffuse )
		{
			diffuse_samples *= diffuse_samples;
			AtSamplerIterator* sampit = AiSamplerIterator(data->diffuse_sampler, sg);
			AiMakeRay(&wi_ray, AI_RAY_DIFFUSE, &sg->P, NULL, AI_BIG, sg);
			AtInt count=0;
			while(AiSamplerGetSample(sampit, samples))
			{
				wi = AiOrenNayarMISSample(dmis, samples[0], samples[1]);
				if (AiV3Dot(wi,sg->Nf) > 0.0f)
				{
					wi_ray.dir = wi;
					AiV3Normalize(H, wi+brdfd.V);
					// attenuate by the fresnel of the top layer
					// we'll assume for now that the specular2 layer does not have a significant attenuation on this
					// layer (it's supposed to be more of a 'mix').
					kt = 1.0f - fresnel(std::max(0.0f,AiV3Dot(H,wi)),eta);
					if (kt > IMPORTANCE_EPS) // only trace a ray if it's going to matter
					{
						AiTrace(&wi_ray, &scrs);
						result_diffuseIndirect +=
						scrs.color*AiOrenNayarMISBRDF(dmis, &wi) / AiOrenNayarMISPDF(dmis, &wi) * kt;
					}
				}
				count++;
			}
			if (count) result_diffuseIndirect /= float(count);
			result_diffuseIndirect *= diffuseColor;
		} // if (do_diffuse)

	} // if (do_diffuse || do_glossy)

	// Emission
	result_emission = emissionColor;

	// Diffusion multiple scattering
	if ( do_sss )
	{
		result_sss = AiSSSPointCloudLookupCubic(sg, sssRadius*sssRadiusColor*sssScale) * diffuseColor;
	}

	// blend sss and direct diffuse
	result_diffuseDirect *= (1-sssMix);
	result_diffuseIndirect *= (1-sssMix);
	result_sss *= sssMix;

	// Refraction
	if (do_transmission)
	{
		result_transmission = beckmannMicrofacetTransmission(sg, sg->N, U, V, wo, data->refraction_sampler,
																transmissionRoughness, transmissionIor,
																sigma_s, sigma_a,
																ssDirection, ssScale, result_ss);
	}

	if (sg->Rt & AI_RAY_CAMERA)
	{
		// write AOVs
		AiAOVSetRGB(sg, "diffuseDirect", result_diffuseDirect);
		AiAOVSetRGB(sg, "multiScatter", result_sss);
		AiAOVSetRGB(sg, "specularDirect", result_glossyDirect);
		AiAOVSetRGB(sg, "specular2Direct", result_glossy2Direct);
		AiAOVSetRGB(sg, "diffuseIndirect", result_diffuseIndirect);
		AiAOVSetRGB(sg, "specularIndirect", result_glossyIndirect);
		AiAOVSetRGB(sg, "specular2Indirect", result_glossy2Indirect);
		AiAOVSetRGB(sg, "singleScatter", result_ss);
		AiAOVSetRGB(sg, "transmission", result_transmission);
		AiAOVSetRGB(sg, "emission", result_emission);

		// write data AOVs
		AtPoint Pref;
		if (!AiUDataGetPnt("Pref", &Pref))
			Pref = sg->Po;
		AtRGB position = AiColorCreate(sg->P.x, sg->P.y, sg->P.z);
		AtRGB referenceposition = AiColorCreate(Pref.x, Pref.y, Pref.z);
		AiAOVSetRGB(sg, "position", position);
		AiAOVSetRGB(sg, "referenceposition", referenceposition);
		AtRGB normal = AiColorCreate(sg->Nf.x, sg->Nf.y, sg->Nf.z);
		AiAOVSetRGB(sg, "normal", normal);
		AtRGB uv = AiColorCreate(sg->u, sg->v, 0.0f);
		AiAOVSetRGB(sg, "uv", uv);
		AtRGB depth = AiColorCreate(sg->Rl, AiV3Dot(sg->Nf, wo), 0.0f);
	}

	// Sum final result from temporaries
	//
	sg->out.RGB =  	 result_diffuseDirect
					+result_sss
					+result_glossyDirect
					+result_glossy2Direct
					+result_diffuseIndirect
					+result_glossyIndirect
					+result_glossy2Indirect
					+result_ss
					+result_transmission
					+result_emission;
}
