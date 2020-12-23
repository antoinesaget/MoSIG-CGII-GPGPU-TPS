#version 410
#define M_PI 3.14159265358979323846
#define EPSILON 0.0001

uniform mat4 mat_inverse;
uniform mat4 persp_inverse;
uniform sampler2D envMap;
uniform vec3 center;
uniform float radius;

uniform bool transparent;
uniform float shininess;
uniform float eta;

in vec4 position;

out vec4 fragColor;

// From https://en.wikipedia.org/wiki/UV_mapping
vec4 getColorFromEnvironment(in vec3 direction)
{
    vec2 uv = vec2(0.5 + atan(direction.z, direction.x)/(2*M_PI), 0.5 + asin(direction.y)/M_PI);
    vec3 color = texture(envMap, uv).rgb;
    return vec4(color, 1);
}

bool raySphereIntersect(in vec3 start, in vec3 direction, in vec3 sphereCenter, in float sphereR, out vec3 newPoint) {
    vec3 hypothenuse = sphereCenter - start;
    float t = dot(direction, hypothenuse);
    float d2 = dot(hypothenuse, hypothenuse) - t*t;
    float r2 = sphereR*sphereR;
    
    if (d2 <= r2) {
        float k = sqrt(r2 - d2);
        if (t - k > EPSILON) { // If intersection at t1
            newPoint = start + direction * (t-k);
        } else if (t + k > EPSILON) { // If intersection at t2
            newPoint = start + direction * (t+k);
        } else { // If intersection behind us
            return false;
        }
        return true;
    } else { // If no intersection
        return false;
    }
}

/*
This function is inspired from the class notes. 
It return a single float and take a single real eta as a parameter.
Here, for the function to work, eta must be eta_out/eta_in
*/
float fresnel_real(vec3 half, vec3 outVect, float eta) {
	float cosTheta = dot(half, outVect);
	float sinTheta = length(cross(half, outVect));

	float ci = sqrt(max(eta*eta - sinTheta*sinTheta, 0.0));
	float Fs = pow((cosTheta - ci) / (cosTheta + ci), 2);
	float Fp = pow((eta*eta*cosTheta - ci) / (eta*eta*cosTheta + ci), 2);
	
     float F = (Fs + Fp) / 2.0;
	return F;
}

/*
Compute the color of a pixel at a given intersection p of the sphere.
If the sphere is transparent then bounces are computed.
If the sphere is opaque we just mix the color of the sphere with  the reflected ray.

No sure of the behaviour of this function for eta < 1. I tried to understand and managed to get a little understanding of what's supposed to happen thanks to https://phet.colorado.edu/sims/html/bending-light/latest/bending-light_fr.html but I think I've done something wrong.
*/
vec4 colorAtIntersection(vec3 p, vec3 In, vec3 center, float radius, float eta1, float eta2) {
    if (!transparent) { // When the sphere is opaque
        vec4 C = vec4(1); // Color of the sphere
        vec3 n = normalize(p-center);
        vec3 reflected = reflect(In, n);
        float F = fresnel_real(n, reflected, eta2/eta1);
        // fcoef = 1.0; // Uncomment to get a classic 100% reflexive metal sphere

        return mix(C, getColorFromEnvironment(reflect(In, n)), F); // The color of the pixel is a mix between the color of the sphere and the reflected ray modulated by the fresnel coefficient.
    } else { // When the sphere is transparent we have to compute multiple bounces
        // TODO
    }
}

void main(void)
{
    // Step 1: I need pixel coordinates. Division by w?
    vec4 worldPos = position;
    worldPos.z = 1; // near clipping plane
    worldPos = persp_inverse * worldPos;
    worldPos /= worldPos.w;
    worldPos.w = 0;
    worldPos = normalize(worldPos);
    // Step 2: ray direction:
    vec3 u = normalize((mat_inverse * worldPos).xyz);
    vec3 eye = (mat_inverse * vec4(0, 0, 0, 1)).xyz;

    vec3 intersection = vec3(0);
    bool isIntersect = raySphereIntersect(eye, u, center, radius, intersection);

    vec4 resultColor = vec4(0);
    if (isIntersect) { // Intersection with the sphere, then compute the color of the pixel at that intersection
        vec3 In = normalize(intersection - eye);
        resultColor = colorAtIntersection(intersection, In, center, radius, 1.0, eta);
    } else { // No interesection, then get the pixel from the environment map
        resultColor = getColorFromEnvironment(u);
    }

    fragColor = resultColor;
}
