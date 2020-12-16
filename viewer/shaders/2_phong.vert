#version 410

uniform mat4 matrix;
uniform mat4 perspective;
uniform mat3 normalMatrix;
uniform bool noColor;
uniform vec3 lightPosition;

// World coordinates
in vec4 vertex;
in vec4 normal;
in vec4 color;

// Camera-space coordinates
out vec4 eyeVector;
out vec4 lightVector;
out vec4 lightSpace; // placeholder for shadow mapping
out vec4 vertColor;
out vec4 vertNormal;

void main( void )
{
    if (noColor) vertColor = vec4(0.2, 0.6, 0.7, 1.0 );
    else vertColor = color;

    vertNormal.xyz = normalize(normalMatrix * normal.xyz);
    vertNormal.w = 0.0;

	vec4 vertex_cameraSpace = matrix * vertex;
	lightVector = matrix * vec4(lightPosition,1.0) - vertex_cameraSpace;
	eyeVector = vec4(0, 0, 0, 1) - vertex_cameraSpace;

	gl_Position = perspective * vertex_cameraSpace;
}
