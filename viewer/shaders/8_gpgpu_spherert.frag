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

    vec4 resultColor = isIntersect ? vec4(1): getColorFromEnvironment(u);
    fragColor = resultColor;
}
