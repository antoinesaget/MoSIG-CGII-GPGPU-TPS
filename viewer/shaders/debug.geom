#version 410
layout (triangles) in;
layout (line_strip, max_vertices = 12) out;

uniform mat4 perspective;
uniform float normalLength;

in vec4 vertNormal[];
in vec4 lightVector[];

out vec4 vertColor;

/*
This geometry can be adapted by the developers for debugging
and will be overlayed on top of the render.
Here for example it's set to create a line for each normal.

This geometry shader was particularly usefull to debug the different vectors (vector to light, normal, half vector, etc.)

The vertex shader input is the currently selected shader by the user, the fragment shader output is debug.frag
*/

void GenerateLine(int index)
{
	vec4 n = normalize(vertNormal[index]);
    vertColor = vec4(1.0, 0.0, 0.0, 1.0);

    gl_Position = gl_in[index].gl_Position;
    EmitVertex();
    
    gl_Position = gl_in[index].gl_Position + normalize(perspective * n) * normalLength;
    EmitVertex();    
    EndPrimitive();

	vec4 L = normalize(lightVector[index]);
    vertColor = vec4(0.0, 0.0, 1.0, 1.0);

    gl_Position = gl_in[index].gl_Position;
    EmitVertex();
    
    gl_Position = gl_in[index].gl_Position + normalize(perspective * L) * normalLength * 5;
    EmitVertex();    
    EndPrimitive();
}

void main()
{
    GenerateLine(0); // first vertex normal
    GenerateLine(1); // second vertex normal
    GenerateLine(2); // third vertex normal
}
