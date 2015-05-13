#include "main.h"

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "control.h"


/*
 * IMPORTANT - DO NOT CHANGE THIS FILE - IMPORTANT
 */


int  window_width  = 300;
int  window_height = 300;

Image* currentImage  = NULL;
Image* originalImage = NULL;
ImageGraph* currentImageGraph = NULL;
vector<vector<bool>>* similarityGraph = NULL;
vector<Curve> * curveVector = NULL;

bool quietMode = false;
bool textMode  = false;


int main (int argc, char** argv)
{
    // initialize parameters
    char* toLoad = init(argc, argv);
    
    if (textMode)
    {
        if (toLoad)
            image_load(toLoad);
        textMenuLoop();
    }
    else
    {
        // set up the window
        glutInit(&argc, &argv[0]);
        glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
        glutInitWindowPosition(100,100);
        glutInitWindowSize(window_width, window_height);
        glutCreateWindow("hmc cs155 image processor");
        
        // register call back functions
        glutDisplayFunc(display);
        glutReshapeFunc(unreshape);
        
        glClearColor(0.0,0.0,0.0,0.0);
        glDisable(GL_DEPTH_TEST);
        
        // setup main menu
        make_menu();
        
        // register keyboard callback function
        glutKeyboardFunc(keyboard_func);
        
        if (toLoad)
            image_load(toLoad);
        
        // wait for something to happen
        glutMainLoop();
    }
    return 0;
}


char* init (int argc, char** argv)
{
    // init random number generator
    //srand48(time(0));
    
    char* toLoad = NULL;
    
    // parse the command line options
    bool noMoreArgs  = false;
    bool noMoreFlags = false;
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            if (noMoreArgs)
                usage();
            
            if (!noMoreFlags && argv[i][0] == '-')
            {
                switch (argv[i][1])
                {
                    case 't':
                        textMode = true;
                        break;
                        
                    case 'q':
                        quietMode = true;
                        break;
                        
                    case '-':
                        if (argv[i][2] == '\0')
                            noMoreFlags = true;
                        else
                            usage();
                        break;
                        
                    default:
                        usage();
                }
            }
            else
            {
                noMoreArgs = true;
                //        image_load(argv[i]);
                toLoad = argv[i];
            }
        }
    }
    
    return toLoad;
}


void usage ()
{
    cerr << "usage: ./ip [ -t ] [ -q ] [ -- ] [ file ]" << endl;
    exit(-1);
}


float N (int i, int k, float t, vector<int>& knots) {
    if (k == 1) {
        if (knots[i] <= t && t < knots[i+1]) return 1;
        else return 0;
        
    }
    return ((float) t - i) / ((float) k) * N(i, k - 1, t, knots) +
    ((float) i + k - t + 1)/((float) k) * N(i+1, k-1, t, knots);
}

void drawCurve(int startPoint, Curve & curve) {
    vector<vector<GLfloat>> points = curve.verticies;
    long numPoints = points.size();

    vector<int> knots(numPoints + 2);
    for (int i = 0; i < numPoints + 2; ++i) {
        knots[i] = i;
    }
    if (points.size() < 2) return;
    points.push_back(points[points.size()-1]);
    /*  approximate the curve by a line strip through sample points	*/
    //glEnable(GL_LINE_WIDTH);
    glLineWidth(3.f);
    glPointSize(5);
    glPushMatrix();
    glColor3f(0,1,0);
    
    glBegin(GL_POINTS);
    for (int i = 0; i < points.size(); ++i) {
        GLfloat poly[2] = { points[i][0], points[i][1]};
        glVertex3f(points[i][0], points[i][1], 0);
    }
    glEnd();
    
    
    GLfloat translate[2] = {points[0][0], points[0][1]};
    glTranslatef(translate[0], translate[1], 0.);
    
    glColor3f(1,0,0);
    
    glBegin(GL_LINE_STRIP);
    float numSamples=100.;
    float t=0;
    
    
    // move to origin, translate later.
    
    for (int i = 0; i < numPoints; ++i)
        for (int k = 0; k < 2; ++k) points[i][k] -= translate[k];
    /*
    
    while (t < numPoints) {
        float polyVal[3] = {0., 0., 0.};
        for (int i = 0; i < numPoints; ++i) {
            for (int k = 0; k < 2; ++k) polyVal[k] += N(i, 3, t, knots) * points[i][k];
        }
        glVertex3fv(polyVal);
        t += ((float) numPoints)/numSamples;
    }*/
    
    
    /* the curve ends at a control point when t=1  				*/
    /* because the increment 1.0/numSamples  has finite precision	*/
    /* t probably won't hit 1.0 exactly, so we force it			*/
    
    glEnd();
    //glDisable(GL_LINE_WIDTH);
    glPopMatrix();
    
    for (int i = 0; i < numPoints; ++i)
        for (int k = 0; k < 2; ++k) points[i][k] += translate[k];
}

void display ()
{
    if (textMode)
        return;
    
    // check if there have been any openGL problems
    GLenum errCode = glGetError();
    if (errCode != GL_NO_ERROR)
    {
        const GLubyte* errString = gluErrorString(errCode);
        cerr << "OpenGL error: " << errString << endl;
    }
    
    // clear the frame buffer
    glClear(GL_COLOR_BUFFER_BIT);
    
    // draw the image
    if (currentImageGraph) {
        for (int i = 0; i < currentImageGraph->size(); ++i)
            currentImageGraph->at(i).glDrawPolygonWrapper();
        if (curveVector) {
            for (int i = 0; i < curveVector->size(); ++i) {
                drawCurve(0, curveVector->at(i));
            }
        }
    }
    else if (currentImage)
        currentImage->glDrawPixelsWrapper();
    
    // swap buffers
    //currentImage->glDrawPixelsWrapper();
    glutSwapBuffers();
}


void unreshape (int width, int height)
{
    // don't allow user to manuall resize the window
    reshape(window_width, window_height);
}


void reshape (int width, int height)
{
    // set window height and width
    window_width  = max(width,  64);
    window_height = max(height, 64); 
    
    if (textMode)
        return;
    
    // change the actual window's size
    glutReshapeWindow(window_width, window_height);
    
    // the lower left corner of the viewport is 0,0
    // the upper right corner is width, height
    glViewport(0, 0, (GLint) window_width, (GLint) window_height);  
    
    // setup orthographic projection matrix
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0.0, window_width, 0.0, window_height);
    
    // default mode should be modelview
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
