/**************************************
*                                     *
*   Jeff Molofee's Basecode Example   *
*          nehe.gamedev.net           *
*                2001                 *
*                                     *
**************************************/

#include <windows.h>												// Header File For Windows
#include <gl\gl.h>													// Header File For The OpenGL32 Library
#include <gl\glu.h>													// Header File For The GLu32 Library
#include "NeHeGL.h"													// Header File For NeHeGL

#include "math.h"												    // NEW: Needed For Sqrtf
#include "ArcBall.h"												// NEW: ArcBall Header
#include "glwindows.h"
#include "TriangleMesh.h"

#pragma comment( lib, "opengl32.lib" )								// Search For OpenGL32.lib While Linking
#pragma comment( lib, "glu32.lib" )									// Search For GLu32.lib While Linking
//#pragma comment( lib, "glaux.lib" )									// Search For GLaux.lib While Linking

#ifndef CDS_FULLSCREEN												// CDS_FULLSCREEN Is Not Defined By Some
#define CDS_FULLSCREEN 4											// Compilers. By Defining It This Way,
#endif																// We Can Avoid Errors

GL_Window*	g_window;
Keys*		g_keys;

// User Defined Variables
GLUquadricObj *quadratic;											// Used For Our Quadric

const float PI2 = 2.0*3.1415926535f;								// PI Squared

Matrix4fT   Transform   = {  1.0f,  0.0f,  0.0f,  0.0f,				// NEW: Final Transform
                             0.0f,  1.0f,  0.0f,  0.0f,
                             0.0f,  0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  0.0f,  1.0f };

Matrix3fT   LastRot     = {  1.0f,  0.0f,  0.0f,					// NEW: Last Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

Matrix3fT   ThisRot     = {  1.0f,  0.0f,  0.0f,					// NEW: This Rotation
                             0.0f,  1.0f,  0.0f,
                             0.0f,  0.0f,  1.0f };

ArcBallT    ArcBall(640.0f, 480.0f);				                // NEW: ArcBall Instance
Point2fT    MousePt;												// NEW: Current Mouse Point
bool        isClicked  = false;										// NEW: Clicking The Mouse?
bool        isRClicked = false;										// NEW: Clicking The Right Mouse Button?
bool        isDragging = false;					                    // NEW: Dragging The Mouse?

bool        isMClicked = false;
float zoomSize = 1.0;
Point2fT    LastMPt;	

float LastMoveX = 0.0f;
float LastMoveY = 0.0f;

float NowMoveX = 0.0f;
float NowMoveY = 0.0f;

float g_width = 640.0f;
float g_height = 480.0f;


BOOL Initialize (GL_Window* window, Keys* keys)						// Any GL Init Code & User Initialiazation Goes Here
{
	g_window	= window;
	g_keys		= keys;

	// Start Of User Initialization
    isClicked   = false;								            // NEW: Clicking The Mouse?
    isDragging  = false;							                // NEW: Dragging The Mouse?

	glClearColor (0.0f, 0.0f, 0.0f, 0.5f);							// Black Background
	glClearDepth (1.0f);											// Depth Buffer Setup
	glDepthFunc (GL_LEQUAL);										// The Type Of Depth Testing (Less Or Equal)
	glEnable (GL_DEPTH_TEST);										// Enable Depth Testing
	glShadeModel (GL_FLAT);											// Select Flat Shading (Nice Definition Of Objects)
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);				// Set Perspective Calculations To Most Accurate

	quadratic=gluNewQuadric();										// Create A Pointer To The Quadric Object
	gluQuadricNormals(quadratic, GLU_SMOOTH);						// Create Smooth Normals
	gluQuadricTexture(quadratic, GL_TRUE);							// Create Texture Coords

	glEnable(GL_LIGHT0);											// Enable Default Light
	glEnable(GL_LIGHTING);											// Enable Lighting

	glEnable(GL_COLOR_MATERIAL);									// Enable Color Material

	return TRUE;													// Return TRUE (Initialization Successful)
}

void Deinitialize (void)											// Any User DeInitialization Goes Here
{
	gluDeleteQuadric(quadratic);
}

void Update (DWORD milliseconds)									// Perform Motion Updates Here
{
	if (g_keys->keyDown [VK_ESCAPE] == TRUE)						// Is ESC Being Pressed?
		TerminateApplication (g_window);							// Terminate The Program

	if (g_keys->keyDown [VK_F1] == TRUE)							// Is F1 Being Pressed?
		ToggleFullscreen (g_window);								// Toggle Fullscreen Mode

    if (isRClicked)													// If Right Mouse Clicked, Reset All Rotations
    {
		Matrix3fSetIdentity(&LastRot);								// Reset Rotation
		Matrix3fSetIdentity(&ThisRot);								// Reset Rotation
        Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);		// Reset Rotation

		LastMoveX = 
		NowMoveX =
        LastMoveY = 
		NowMoveY = 0.0f;
    }

    if (!isDragging)												// Not Dragging
    {
        if (isClicked)												// First Click
        {
			isDragging = true;										// Prepare For Dragging
			LastRot = ThisRot;										// Set Last Static Rotation To Last Dynamic One
			ArcBall.click(&MousePt);								// Update Start Vector And Prepare For Dragging

			LastMoveX = NowMoveX;
			LastMoveY = NowMoveY;
			LastMPt = MousePt;
        }
    }
    else
    {
        if (isClicked && !isMClicked)												// Still Clicked, So Still Dragging
        {
            Quat4fT     ThisQuat;

            ArcBall.drag(&MousePt, &ThisQuat);						// Update End Vector And Get Rotation As Quaternion
            Matrix3fSetRotationFromQuat4f(&ThisRot, &ThisQuat);		// Convert Quaternion Into Matrix3fT
            Matrix3fMulMatrix3f(&ThisRot, &LastRot);				// Accumulate Last Rotation Into This One
            Matrix4fSetRotationFromMatrix3f(&Transform, &ThisRot);	// Set Our Final Transform's Rotation From This One
        }else if (isClicked && isMClicked){
			NowMoveX = LastMoveX + (MousePt.s.X - LastMPt.s.X)/g_width;
			NowMoveY = LastMoveY + (MousePt.s.Y - LastMPt.s.Y)/g_height;
		}
        else														// No Longer Dragging
            isDragging = false;
    }
}

void Torus(float MinorRadius, float MajorRadius)					// Draw A Torus With Normals
{
	int i, j;
	glBegin( GL_TRIANGLE_STRIP );									// Start A Triangle Strip
		for (i=0; i<20; i++ )										// Stacks
		{
			for (j=-1; j<20; j++)									// Slices
			{
				float wrapFrac = (j%20)/(float)20;
				float phi = PI2*wrapFrac;
				float sinphi = float(sin(phi));
				float cosphi = float(cos(phi));

				float r = MajorRadius + MinorRadius*cosphi;

				glNormal3f(float(sin(PI2*(i%20+wrapFrac)/(float)20))*cosphi, sinphi, float(cos(PI2*(i%20+wrapFrac)/(float)20))*cosphi);
				glVertex3f(float(sin(PI2*(i%20+wrapFrac)/(float)20))*r,MinorRadius*sinphi,float(cos(PI2*(i%20+wrapFrac)/(float)20))*r);

				glNormal3f(float(sin(PI2*(i+1%20+wrapFrac)/(float)20))*cosphi, sinphi, float(cos(PI2*(i+1%20+wrapFrac)/(float)20))*cosphi);
				glVertex3f(float(sin(PI2*(i+1%20+wrapFrac)/(float)20))*r,MinorRadius*sinphi,float(cos(PI2*(i+1%20+wrapFrac)/(float)20))*r);
			}
		}
	glEnd();														// Done Torus
}

void Draw (Rgbd::GLMesh* mesh)
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				// Clear Screen And Depth Buffer
	glLoadIdentity();												// Reset The Current Modelview Matrix
	glTranslatef(-1.5f,0.0f,-6.0f);									// Move Left 1.5 Units And Into The Screen 6.0

	glScalef(zoomSize,zoomSize,zoomSize);
	glTranslatef (NowMoveX*5, -NowMoveY*5, 0.0);
    glPushMatrix();													// NEW: Prepare Dynamic Transform
    glMultMatrixf(Transform.M);										// NEW: Apply Dynamic Transform
	
	if(mesh->gl_meshes.size() > 0){
		glBegin(GL_TRIANGLES);
		for(std::map<std::tuple<int,int,int>,Rgbd::GLTriangle>::iterator it = mesh->gl_meshes.begin();it != mesh->gl_meshes.end();it++){
			//glBegin(GL_LINE_LOOP);
			Rgbd::GLTriangle m = it->second;
			Rgbd::GLVertex v[3];
			v[0] = mesh->gl_vertexes.at(m.vertexes[0]);
			v[1] = mesh->gl_vertexes.at(m.vertexes[1]);
			v[2] = mesh->gl_vertexes.at(m.vertexes[2]);
			glNormal3f(v[0].v[3], v[0].v[4], v[0].v[5]);
			glVertex3f(v[0].v[0], v[0].v[1], v[0].v[2]);
			glNormal3f(v[1].v[3], v[1].v[4], v[1].v[5]);
			glVertex3f(v[1].v[0], v[1].v[1], v[1].v[2]);
			glNormal3f(v[2].v[3], v[2].v[4], v[2].v[5]);
			glVertex3f(v[2].v[0], v[2].v[1], v[2].v[2]);
			//glEnd();
		}
		glEnd();
	}
	//
	//else if(Sdf::sdf != NULL){
	//	glBegin(GL_POINTS);
	//	glColor3f(1.0,1.0,1.0); //设置点颜色
	//	Sdf::sdf->TraversalDraw();
	//	glEnd();
	//	
	//}else if(Sdf::tsdf != NULL){
	//	glBegin(GL_POINTS);
	//	glColor3f(1.0,1.0,1.0); //设置点颜色
	//	Sdf::tsdf->TraversalDraw();
	//	glEnd();
	//	
	//}
	

	

	/*
	if(Sdf::sdf->gen != NULL){
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < Sdf::sdf->gen->m_nTriangles; i++) {
			VECTOR3D vec1, vec2, normal;
			unsigned int id[3];
			id[0] = Sdf::sdf->gen->m_piTriangleIndices[i*3];
			id[1] = Sdf::sdf->gen->m_piTriangleIndices[i*3+1];
			id[2] =Sdf::sdf->gen->m_piTriangleIndices[i*3+2];
			for(int j = 0;j < 3;j++){
				glNormal3f(Sdf::sdf->gen->m_pvec3dNormals[id[j]][0],Sdf::sdf->gen->m_pvec3dNormals[id[j]][1],Sdf::sdf->gen->m_pvec3dNormals[id[j]][2]);
				glVertex3f(Sdf::sdf->gen->m_ppt3dVertices[id[j]][0],Sdf::sdf->gen->m_ppt3dVertices[id[j]][1],Sdf::sdf->gen->m_ppt3dVertices[id[j]][2]);
			}
			
		}
		glEnd();
	}
	*/
	
	
	
	
    glPopMatrix();													// NEW: Unapply Dynamic Transform

	glFlush ();														// Flush The GL Rendering Pipeline
}


/*
int DrawGLScene(GLvoid)									// Here's Where We Do All The Drawing
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);	// Clear Screen And Depth Buffer
	glLoadIdentity();									// Reset The Current Modelview Matrix
	glTranslatef(-1.5f,0.0f,-6.0f);						// Move Left 1.5 Units And Into The Screen 6.0
	
	glRotatef(rtri,0.0f,1.0f,0.0f);
	glBegin(GL_TRIANGLES);							// Start Drawing A Triangle
	for(std::vector<GLMesh>::iterator it = GL_Mesh::gl_meshes->begin();it != GL_Mesh::gl_meshes->end();it++){
		GLMesh m = *it;
		glNormal3f(m.mesh[0][3], m.mesh[0][4], m.mesh[0][5]);
		glVertex3f(m.mesh[0][0], m.mesh[0][1], m.mesh[0][2]);
		glNormal3f(m.mesh[1][3], m.mesh[1][4], m.mesh[1][5]);
		glVertex3f(m.mesh[1][0], m.mesh[1][1], m.mesh[1][2]);
		glNormal3f(m.mesh[2][3], m.mesh[2][4], m.mesh[2][5]);
		glVertex3f(m.mesh[2][0], m.mesh[2][1], m.mesh[2][2]);
	}
    glEnd();	
	rtri+=3.0f;	

	glEnable(GL_LIGHTING);
	return TRUE;										// Keep Going
}
*/
