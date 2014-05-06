//
//  main.cpp
//  SPHSim_C++
//
//  Created by Peihong Guo on 4/22/14.
//  Copyright (c) 2014 Peihong Guo. All rights reserved.
//

#include <iostream>
#include "SPHSystem.h"
#include "glfw3.h"
#include "OpenGL/glu.h"
#include "glm.hpp"
#include "gtc/matrix_transform.hpp"


static void error_callback(int error, const char* description)
{
    fputs(description, stderr);
}
static void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(window, GL_TRUE);
}

void drawCube() {
    // Render a cube
    glBegin( GL_LINE_LOOP );
    // Top face
    glColor3f(   0.0f, 1.0f,  0.0f );  // Green
    glVertex3f(  1.0f, 1.0f, -1.0f );  // Top-right of top face
    glVertex3f( -1.0f, 1.0f, -1.0f );  // Top-left of top face
    glVertex3f( -1.0f, 1.0f,  1.0f );  // Bottom-left of top face
    glVertex3f(  1.0f, 1.0f,  1.0f );  // Bottom-right of top face
    glEnd();
    
    glBegin( GL_LINE_LOOP );
    // Bottom face
    glColor3f(   1.0f,  0.5f,  0.0f ); // Orange
    glVertex3f(  1.0f, -1.0f, -1.0f ); // Top-right of bottom face
    glVertex3f( -1.0f, -1.0f, -1.0f ); // Top-left of bottom face
    glVertex3f( -1.0f, -1.0f,  1.0f ); // Bottom-left of bottom face
    glVertex3f(  1.0f, -1.0f,  1.0f ); // Bottom-right of bottom face
    glEnd();
    
    glBegin( GL_LINE_LOOP );
    // Front face
    glColor3f(   1.0f,  0.0f, 0.0f );  // Red
    glVertex3f(  1.0f,  1.0f, 1.0f );  // Top-Right of front face
    glVertex3f( -1.0f,  1.0f, 1.0f );  // Top-left of front face
    glVertex3f( -1.0f, -1.0f, 1.0f );  // Bottom-left of front face
    glVertex3f(  1.0f, -1.0f, 1.0f );  // Bottom-right of front face
    glEnd();
    
    glBegin( GL_LINE_LOOP );
    // Back face
    glColor3f(   1.0f,  1.0f,  0.0f ); // Yellow
    glVertex3f(  1.0f, -1.0f, -1.0f ); // Bottom-Left of back face
    glVertex3f( -1.0f, -1.0f, -1.0f ); // Bottom-Right of back face
    glVertex3f( -1.0f,  1.0f, -1.0f ); // Top-Right of back face
    glVertex3f(  1.0f,  1.0f, -1.0f ); // Top-Left of back face
    glEnd();
    
    glBegin( GL_LINE_LOOP );
    // Left face
    glColor3f(   0.0f,  0.0f,  1.0f);  // Blue
    glVertex3f( -1.0f,  1.0f,  1.0f);  // Top-Right of left face
    glVertex3f( -1.0f,  1.0f, -1.0f);  // Top-Left of left face
    glVertex3f( -1.0f, -1.0f, -1.0f);  // Bottom-Left of left face
    glVertex3f( -1.0f, -1.0f,  1.0f);  // Bottom-Right of left face
    glEnd();
    
    glBegin( GL_LINE_LOOP );
    // Right face
    glColor3f(   1.0f,  0.0f,  1.0f);  // Violet
    glVertex3f(  1.0f,  1.0f,  1.0f);  // Top-Right of left face
    glVertex3f(  1.0f,  1.0f, -1.0f);  // Top-Left of left face
    glVertex3f(  1.0f, -1.0f, -1.0f);  // Bottom-Left of left face
    glVertex3f(  1.0f, -1.0f,  1.0f);  // Bottom-Right of left face
    glEnd();
}

void drawParticles(SPHSystem* sph) {
    auto p = sph->particles();
    glPointSize(2.0f);
    glColor3f(1, 1, 1);
    float r, g, b;
    glBegin(GL_POINTS);
    for(auto x : p) {
        auto pos = x*2.0f-1.0f;
        r = x.y/0.25, b = 1.0 - r, g = 0.0;
        glColor3f(r, g, b);
        glVertex3f(pos.x, pos.y, pos.z);
    }
    glEnd();
}

int main(int argc, const char * argv[])
{
    SPHSystem *sph;
    GLFWwindow* window;
    
    if( argc > 1 )
        sph = new SPHSystem(argv[1]);
    else
        sph = new SPHSystem;
    
    sph->init();
    
    glfwSetErrorCallback(error_callback);
    if (!glfwInit())
        exit(EXIT_FAILURE);
    window = glfwCreateWindow(640, 480, "SPHSim", NULL, NULL);
    if (!window)
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }
    glfwMakeContextCurrent(window);
    glfwSetKeyCallback(window, key_callback);
    while (!glfwWindowShouldClose(window))
    {
        sph->step();
        
        // setup viewing params
        float ratio;
        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        ratio = width / (float) height;
        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
        
        glm::mat4 mproj = glm::perspective(45.0f, ratio, 0.0001f, 10.0f);
        glm::mat4 mview = glm::lookAt(glm::vec3(0, 0.0, -3.5), glm::vec3(0), glm::vec3(0, 1, 0));
        glm::mat4 mmodel = glm::mat4(1.0f);
        glm::mat4 mmv = mview * mmodel;
        
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glMultMatrixf(&mproj[0][0]);
        
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        glMultMatrixf(&mmv[0][0]);
        
        glShadeModel(GL_SMOOTH);
        glEnable( GL_DEPTH_TEST );
        glDepthFunc(GL_LESS);

        // render the scene
        drawParticles(sph);
        drawCube();
        
        glfwSwapBuffers(window);
        glfwPollEvents();
    }
    glfwDestroyWindow(window);
    glfwTerminate();
    exit(EXIT_SUCCESS);
    
    std::cout << "Done!\n";
    return 0;
}

