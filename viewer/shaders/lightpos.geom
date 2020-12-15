#version 410
layout (points) in;
layout (line_strip, max_vertices = 2) out;

uniform vec3 lightPosition;
uniform float ground_level;

uniform mat4 matrix;
uniform mat4 perspective;

/*
Simple geometry shader that create a line from the light source to the ground
for debug purposes.
*/ 

void main()
{
    gl_Position = perspective * matrix * vec4(lightPosition, 1.0);
    EmitVertex();
    
    gl_Position = perspective * matrix * vec4(lightPosition.x, ground_level, lightPosition.z, 1.0);
    EmitVertex();
    
    EndPrimitive();
}
