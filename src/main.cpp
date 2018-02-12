#include "gl.h"
#include <GLFW/glfw3.h>

#include <cmath>
#include <vector>
#include <cassert>
#include <iostream>
#include <sstream>
#include <map>
#include <algorithm>

#include <vecmath.h>
#include "camera.h"
#include "starter0_util.h"
#include "teapot.h"
#include "tuple.h"
#include "vertexrecorder.h"

using namespace std;

// Globals
//uint32_t program;
Camera camera;
bool gMousePressed = false;
int subsurfed = 0;
bool wireframe = false;
GLuint program;
GLuint program_color;

// Variables beginning with "ss" represent subsurfed equivalents of the variable
// This is the list of points (3D vectors)
vector<Vector3f> vecv;
vector<Vector3f> ssvecv;

// This is the list of normals (also 3D vectors)
vector<Vector3f> vecn;
vector<Vector3f> ssvecn;

// This is the list of faces (indices into vecv and vecn)
vector<vector<int>> vecf;
vector<vector<int>> ssvecf;

// You will need more global variables to implement color and position changes
int colorCount = 0;
float lightX = 2.0;
float lightY = 3.0;
const float pi = atanf(1) * 4;
const float increment = 45 * pi / 180;
float currentAngle = 0;

void keyCallback(GLFWwindow* window, int key,
    int scancode, int action, int mods)
{
    if (action == GLFW_RELEASE) { // only handle PRESS and REPEAT
        return;
    }

    // Special keys (arrows, CTRL, ...) are documented
    // here: http://www.glfw.org/docs/latest/group__keys.html
    if (key == GLFW_KEY_ESCAPE) {
        glfwSetWindowShouldClose(window, GLFW_TRUE);
    } else if (key == 'W') {
        wireframe = !wireframe;
	} else if (key == GLFW_KEY_C) {
		colorCount++;
		colorCount = colorCount % 4; // Just in case the counter exceeds the length of the color array
	} else if (key == GLFW_KEY_UP) {
		lightY += 0.5;
	} else if (key == GLFW_KEY_DOWN) {
		lightY -= 0.5;
	} else if (key == GLFW_KEY_LEFT) {
		lightX -= 0.5;
	} else if (key == GLFW_KEY_RIGHT) {
		lightX += 0.5;
	} else if (key == GLFW_KEY_R) {
		currentAngle += increment;
	} else if (key == GLFW_KEY_S) {
		subsurfed++;
	} else if (key == GLFW_KEY_D) {
		subsurfed= 0;
	} else if (key == GLFW_KEY_SPACE) {
		Matrix4f eye = Matrix4f::identity();
		camera.SetRotation(eye);
		camera.SetCenter(Vector3f(0, 0, 0));
	} else {
        printf("Unhandled key press %d\n", key);
    }
}

static void mouseCallback(GLFWwindow* window, int button, int action, int mods) {
	double xd, yd;
	glfwGetCursorPos(window, &xd, &yd);
	int x = (int)xd;
	int y = (int)yd;

	int lstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT);
	int rstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT);
	int mstate = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE);
	if (lstate == GLFW_PRESS) {
		gMousePressed = true;
		camera.MouseClick(Camera::LEFT, x, y);
	} else if (rstate == GLFW_PRESS) {
		gMousePressed = true;
		camera.MouseClick(Camera::RIGHT, x, y);
	} else if (mstate == GLFW_PRESS) {
		gMousePressed = true;
		camera.MouseClick(Camera::MIDDLE, x, y);
	} else {
		gMousePressed = true;
		camera.MouseRelease(x, y);
		gMousePressed = false;
	}
}

static void motionCallback(GLFWwindow* window, double x, double y) {
	if (!gMousePressed) {
		return;
	}
	camera.MouseDrag((int)x, (int)y);
}

void setViewport(GLFWwindow* window) {
	int w, h;
	glfwGetFramebufferSize(window, &w, &h);

	camera.SetDimensions(w, h);
	camera.SetViewport(0, 0, w, h);
	camera.ApplyViewport();
}

void drawTriangle()
{
    // set a reasonable upper limit for the buffer size
    VertexRecorder rec;
    rec.record(Vector3f(0.0, 0.0, 0.0), // Position
        Vector3f(0.0, 0.0, 1.0));// Normal

    rec.record(Vector3f(3.0, 0.0, 0.0),
        Vector3f(0.0, 0.0, 1.0));

    rec.record(Vector3f(3.0, 3.0, 0.0),
        Vector3f(0.0, 0.0, 1.0));
    rec.draw();
}

void drawTeapot()
{
    // set the required buffer size exactly.
    VertexRecorder rec;
    for (int idx : teapot_indices) {
        Vector3f position(teapot_positions[idx * 3 + 0],
            teapot_positions[idx * 3 + 1],
            teapot_positions[idx * 3 + 2]);

        Vector3f normal(teapot_normals[idx * 3 + 0],
            teapot_normals[idx * 3 + 1],
            teapot_normals[idx * 3 + 2]);

        rec.record(position, normal);
    }
    rec.draw();
}

void drawObjMesh() {
    // draw obj mesh here
    // read vertices and face indices from vecv, vecn, vecf
	VertexRecorder rec;
	for (vector<int> face : vecf) {
		if ((int)face.size() == 6) {
			// Read through vecv and vecn using the indices from vecf
			// Indices are in the order [a, c, d, f, g, i] as shown in the instruction document
			rec.record(vecv[face[0] - 1], vecn[face[1] - 1]);
			rec.record(vecv[face[2] - 1], vecn[face[3] - 1]);
			rec.record(vecv[face[4] - 1], vecn[face[5] - 1]);
		} else {
			rec.record(vecv[face[0] - 1], vecn[face[0] - 1]);
			rec.record(vecv[face[1] - 1], vecn[face[1] - 1]);
			rec.record(vecv[face[2] - 1], vecn[face[2] - 1]);
			rec.record(vecv[face[2] - 1], vecn[face[2] - 1]);
			rec.record(vecv[face[3] - 1], vecn[face[3] - 1]);
			rec.record(vecv[face[0] - 1], vecn[face[0] - 1]);
		}
	}
	rec.draw();
}

void ssCalculateNs() {
	// First calculate face normals
	vector<Vector3f> faceNs;
	for (vector<int> face : ssvecf) {
		Vector3f side1 = ssvecv[face[1] - 1] - ssvecv[face[0] - 1];
		Vector3f side2 = ssvecv[face[3] - 1] - ssvecv[face[0] - 1];
		Vector3f norm = Vector3f::cross(side1, side2).normalized();
		faceNs.push_back(norm);
	}
	vector<Vector3f>().swap(ssvecn);
	for (int i = 0; i < (int)ssvecv.size(); i++) {
		Vector3f vertNorm = Vector3f(0.0f);
		for (int f = 0; f < (int)ssvecf.size(); f++) {
			if (find(ssvecf[f].begin(), ssvecf[f].end(), i + 1) != ssvecf[f].end()) {
				vertNorm += faceNs[f];
			}
		}
		vertNorm.normalize();
		ssvecn.push_back(vertNorm);
	}
}

void quadCalculateNs() {
	// First calculate face normals
	vector<Vector3f> faceNs;
	for (vector<int> face : vecf) {
		Vector3f side1 = vecv[face[1] - 1] - vecv[face[0] - 1];
		Vector3f side2 = vecv[face[3] - 1] - vecv[face[0] - 1];
		Vector3f norm = Vector3f::cross(side1, side2).normalized();
		faceNs.push_back(norm);
	}
	vector<Vector3f>().swap(vecn);
	for (int i = 0; i < (int)vecv.size(); i++) {
		Vector3f vertNorm = Vector3f(0.0f);
		for (int f = 0; f < (int)vecf.size(); f++) {
			if (find(vecf[f].begin(), vecf[f].end(), i + 1) != vecf[f].end()) {
				vertNorm += faceNs[f];
			}
		}
		vertNorm.normalize();
		vecn.push_back(vertNorm);
	}
}

void drawSubsurf() {
	// Read vertices from ssvecv, face indices from ssvecf, and calculate normals
	ssCalculateNs();

	VertexRecorder rec;
	for (vector<int> face : ssvecf) {
		rec.record(ssvecv[face[0] - 1], ssvecn[face[0] - 1]);
		rec.record(ssvecv[face[1] - 1], ssvecn[face[1] - 1]);
		rec.record(ssvecv[face[2] - 1], ssvecn[face[2] - 1]);
		rec.record(ssvecv[face[2] - 1], ssvecn[face[2] - 1]);
		rec.record(ssvecv[face[3] - 1], ssvecn[face[3] - 1]);
		rec.record(ssvecv[face[0] - 1], ssvecn[face[0] - 1]);
	}
	rec.draw();

}



// This function is responsible for displaying the object.
void drawScene()
{
	if (subsurfed >= 1) {
		drawSubsurf();
	} else {
		drawObjMesh();
	}
    // drawTeapot();
}

/*void setViewport(GLFWwindow* window)
{
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    // make sure the viewport is square-shaped.
    if (width > height) {
        int offsetx = (width - height) / 2;
        glViewport(offsetx, 0, height, height);
    } else {
        int offsety = (height - width) / 2;
        glViewport(0, offsety, width, width);
    }
}*/

void updateCameraUniforms(GLuint program)
{
    // Set up a perspective view, with square aspect ratio
    float fovy_radians = deg2rad(50.0f);
    float nearz = 1.0f;
    float farz = 100.0f;
    float aspect = 1.0f;
    Matrix4f P = Matrix4f::perspectiveProjection(
        fovy_radians, aspect, nearz, farz);

    // See https://www.opengl.org/sdk/docs/man/html/glUniform.xhtml
    // for the many version of glUniformXYZ()
    // Returns -1 if uniform not found.
    int loc = glGetUniformLocation(program, "P");
    glUniformMatrix4fv(loc, 1, false, P);

    Vector3f eye(0.0, 0.0, 7.0f);
    Vector3f center(0.0, 0.0, 0.0);
    Vector3f up(0.0, 1.0f, -0.2f);
    Matrix4f V = Matrix4f::lookAt(eye, center, up);
    loc = glGetUniformLocation(program, "V");
    glUniformMatrix4fv(loc, 1, false, V);
    loc = glGetUniformLocation(program, "camPos");
    glUniform3fv(loc, 1, eye);

    // Make sure the model is centered in the viewport
    // We translate the model using the "Model" matrix
    Matrix4f M = Matrix4f::translation(0, -2.0, 0);
    loc = glGetUniformLocation(program, "M");
    glUniformMatrix4fv(loc, 1, false, M);

    // Transformation matrices act differently
    // on vectors than on points.
    // The inverse-transpose is what we want.
    Matrix4f N = M.inverse().transposed();
    loc = glGetUniformLocation(program, "N");
    glUniformMatrix4fv(loc, 1, false, N);
}

void updateMaterialUniforms(GLuint program)
{
    // Here are some colors you might use - feel free to add more
    GLfloat diffColors[4][4] = { 
    { 0.5f, 0.5f, 0.9f, 1.0f },
    { 0.9f, 0.5f, 0.5f, 1.0f },
    { 0.5f, 0.9f, 0.3f, 1.0f },
    { 0.3f, 0.8f, 0.9f, 1.0f } };

    // Here we use the first color entry as the diffuse color
    int loc = glGetUniformLocation(program, "diffColor");
    glUniform4fv(loc, 1, diffColors[colorCount]); // Pick from the table

    // Define specular color and shininess
    GLfloat specColor[] = { 0.2f, 0.2f, 0.2f, 1.0f };
    GLfloat shininess[] = { 10.0f };

    // Note that the specular color and shininess can stay constant
    loc = glGetUniformLocation(program, "specColor");
    glUniform4fv(loc, 1, specColor);
    loc = glGetUniformLocation(program, "shininess");
    glUniform1f(loc, shininess[0]);
}

void updateLightUniforms(GLuint program)
{
    // Light Position
    GLfloat lightPos[] = { lightX, lightY, 5.0f, 1.0f };
    int loc = glGetUniformLocation(program, "lightPos");
    glUniform4fv(loc, 1, lightPos);

    // Light Color
    GLfloat lightDiff[] = { 120.0, 120.0, 120.0, 1.0 };
    loc = glGetUniformLocation(program, "lightDiff");
    glUniform4fv(loc, 1, lightDiff);
}

void updateRotationUniforms(GLuint program) {
	GLint model = glGetUniformLocation(program, "M");
	Matrix4f rotation = Matrix4f::rotateY(currentAngle);
	rotation = Matrix4f::translation(0, -2.0, 0) * rotation; // translate in the same way that camera does
	rotation = rotation * camera.GetRotation();
	glUniformMatrix4fv(model, 1, false, rotation);
}

// Function that takes in a string of the form "x/y/z" and returns the array [x, z] where x and z are integers
vector<int> indexStripper(string indices) {
	vector<int> answer;
	string slash = "/";
	// Find each / character and use them to make substrings
	string triple = indices;
	int first = triple.find(slash);
	string firstPiece = triple.substr(0, first);
	answer.push_back(stoi(firstPiece));
	int second = triple.find(slash, first + 1);
	string secondPiece = triple.substr(second + 1);
	answer.push_back(stoi(secondPiece));

	return answer;
}

// If the obj file doesn't specify vector normals, calculate them given the faces
void calculateNs() {
	// First calculate all the face normals
	vector<Vector3f> fn;
	for (vector<int> face : vecf) {
		Vector3f v1 = vecv[face[2] - 1] - vecv[face[0] - 1];
		Vector3f v2 = vecv[face[4] - 1] - vecv[face[0] - 1];
		fn.push_back(Vector3f::cross(v1, v2).normalized());
	}
	// Then for each vertex find which faces it's in and add their normals, normalizing for the final result
	for (int vert = 0; vert < (int) vecv.size(); vert++) {
		Vector3f norm = Vector3f(0.0f);
		vector<int> faces;
		for (int f = 0; f < (int) vecf.size(); f++) {
			if (vecf[f][0] - 1 == vert || vecf[f][2] - 1 == vert || vecf[f][4] - 1 == vert) {
				norm += fn[f];
				faces.push_back(f);
			}
		}
		norm.normalize();
		vecn.push_back(norm);
		for (int face : faces) {
			if (vecf[face][0] - 1 == vert) {
				vecf[face][1] = vecn.size();
			} else if (vecf[face][2] - 1 == vert) {
				vecf[face][3] = vecn.size();
			} else {
				vecf[face][5] = vecn.size();
			}
		}
	}
}

void loadInput()
{
    // load the OBJ file here
	const int MAX_BUFFER_SIZE = 4096;
	char buffer[MAX_BUFFER_SIZE];
	bool calculateNormals = false;
	while (cin.getline(buffer, MAX_BUFFER_SIZE)) {
		stringstream ss(buffer);
		Vector3f v;
		vector<string> vs (4); // For the faces
		vector<int> tri; // Again for the faces
		string s;
		ss >> s;
		//First, read in the vertices, lines start with "v"
		if (s == "v") {
			ss >> v[0] >> v[1] >> v[2]; 
			vecv.push_back(v);
		} else if (s == "vn") { //Now for the normal lines starting with "vn"
			ss >> v[0] >> v[1] >> v[2];
			vecn.push_back(v);
		} else if (s == "f") { // And finally the face lines starting with "f"
			ss >> vs[0] >> vs[1] >> vs[2] >> vs[3];
			char test = vs[3][0];
			// Check if the faces specify normals or not, ie if they have backslashes
			if (vs[0].find('/') != string::npos) {
				// Each piece is of the format "x/y/z", and we only need x and z
				for (int i = 0; i < 3; i++) {
					vector<int> strippedPiece = indexStripper(vs[i]);
					tri.push_back(strippedPiece[0]);
					tri.push_back(strippedPiece[1]);
				}
				vecf.push_back(tri);
			} else if (false) {
				// Each vertex is simply listed, but we have to calculate normals
				tri = { stoi(vs[0]), 0, stoi(vs[1]), 0, stoi(vs[2]), 0 };
				vecf.push_back(tri);
				calculateNormals = true;
			} else {
				tri = { stoi(vs[0]), stoi(vs[1]), stoi(vs[2]), stoi(vs[3]) };
				vecf.push_back(tri);
				quadCalculateNs();
			}
		}
	}
	if (calculateNormals) {
		calculateNs();
	}
}



// Helper function for subsurf that returns a vector of vectors such that vertex i
// belongs to all the faces in result[i]
vector<vector<int>> vertFaces(vector<Vector3f> verts, vector<vector<int>> faces) {
	vector<vector<int>> relation;
	for (int i = 0; i < (int)verts.size(); i++) {
		vector<int> facesBelongedTo;
		for (int j = 0; j < (int)faces.size(); j++) {
			for (int v = 0; v < (int)faces[j].size(); v ++) {
				if (faces[j][v] - 1 == i) {
					facesBelongedTo.push_back(j);
				}
			}
		}
		relation.push_back(facesBelongedTo);
	}
	return relation;
}

// Helper function that iterates through faces and returns a vector of all pairs of vertices
// that have an edge between them
vector<stuple<int, 2>> findEdges(vector<vector<int>> faces) {
	vector<stuple<int, 2>> relation;
	for (vector<int> face : faces) {
		for (int i = 0; i < (int)face.size() - 1; i++) {
			stuple<int, 2> edge = stuple<int, 2>(min(face[i], face[i + 1]), max(face[i], face[i + 1]));
			bool toAdd = true;
			for (int j = 0; j < (int) relation.size(); j++){
				if (relation[j] == edge) {
					toAdd = false;
				}
			}
			if (toAdd) {
				relation.push_back(edge);
			}
		}
		stuple<int, 2> lastEdge = stuple<int, 2>(min(face[(int)face.size() - 1], face[0]), max(face[(int)face.size() - 1], face[0]));
		bool toAdd2 = true;
		for (int k = 0; k < (int)relation.size(); k++) {
			if (relation[k] == lastEdge) {
				toAdd2 = false;
			}
		}
		if (toAdd2) {
			relation.push_back(lastEdge);
		}
	}
	return relation;
}

// Helper function that iterates through the list of edges and faces and matches edges to the faces
// that it borders
map<stuple<int, 2>, stuple<int, 2>> edgeFaces(vector<stuple<int, 2>> edges, vector<vector<int>> vertFaces) {
	map<stuple<int, 2>, stuple<int, 2>> relation;
	for (stuple<int, 2> edge : edges) {
		//cerr << "Edge: <" << edge[0] << ", " << edge[1] << "> to faces: ";
		vector<int> facesToAdd;
		int v1 = edge[0] - 1;
		int v2 = edge[1] - 1;
		vector<int> faces1 = vertFaces[v1];
		sort(faces1.begin(), faces1.end());
		vector<int> faces2 = vertFaces[v2];
		sort(faces2.begin(), faces2.end());
		set_intersection(faces1.begin(), faces1.end(), faces2.begin(), faces2.end(), back_inserter(facesToAdd));
		relation[edge] = stuple<int, 2>(facesToAdd[0], facesToAdd[1]);
		//cerr << "(" << relation[edge][0] << ", " << relation[edge][1] << ")" << endl;
	}
	return relation;
}

// Helper function that calculates all the edge points
map<stuple<int, 2>, Vector3f> edgePoints(vector<stuple<int, 2>> edges, vector<Vector3f> verts,
	vector<Vector3f> facePoints, map<stuple<int, 2>, stuple<int, 2>> edgeFaces) {
	map<stuple<int, 2>, Vector3f> relation;
	for (stuple<int, 2> edge : edges) {
		Vector3f v1 = verts[edge[0] - 1];
		Vector3f v2 = verts[edge[1] - 1];
		stuple<int, 2> faces = edgeFaces[edge];
		Vector3f v3 = facePoints[faces[0]];
		Vector3f v4 = facePoints[faces[1]];
		Vector3f ave = v1 + v2 + v3 + v4;
		ave /= 4;
		relation[edge] = ave;
	}
	return relation;
}

// Calculates the midpoints of edges and adds them to a map
map<stuple<int, 2>, Vector3f> edgeMids(vector<Vector3f> verts, vector<stuple<int, 2>> edges) {
	map<stuple<int, 2>, Vector3f> relation;
	for (stuple<int, 2> edge : edges) {
		Vector3f v1 = verts[edge[0] - 1];
		Vector3f v2 = verts[edge[1] - 1];
		Vector3f mid = v1 + v2;
		mid /= 2;
		relation[edge] = mid;
	}
	return relation;
}


// Helper function that calculates the average of the midpoints of edges that given vertices belong to
vector<Vector3f> aveMidEdges(vector<Vector3f> verts, vector<stuple<int, 2>> edges) {
	vector<Vector3f> averages;
	map<stuple<int, 2>, Vector3f> midpts = edgeMids(verts, edges);
	for (int i = 0; i < (int)verts.size(); i++) {
		Vector3f aves = Vector3f(0.0f);
		int divisor = 0;
		for (int j = 0; j < (int)edges.size(); j++) {
			if (edges[j][0] - 1 == i || edges[j][1] - 1 == i) {
				aves += midpts[edges[j]];
				divisor++;
			}
		}
		aves /= (float)divisor;
		averages.push_back(aves);
	}
	return averages;
}

// Applies the Catmull-Clark subdivision surface algorithm
void subsurfCC(vector<Vector3f> verts, vector<vector<int>> faces, vector<Vector3f> &newPoints, vector<vector<int>> &newFaces, int reps = 1) {
	vector<Vector3f> facePoints;
	vector<vector<int>> vertFs = vertFaces(verts, faces); // Vertex i belongs to the faces in vertFaces[i]
	vector<stuple<int, 2>> edges = findEdges(faces);
	map<stuple<int, 2>, stuple<int, 2>> edgeFs = edgeFaces(edges, vertFs); // The keys are the edges (represented by
	// the two vertices that form the endpoints, and the values are the two faces that the 
	// edge borders
	// First, create face points by averaging the vertices of every face, and calculate face normals
	for (vector<int> face : faces) {
		int divisor = (int)face.size();
		Vector3f ave = Vector3f(0.0f);
		for (int i = 0; i < (int)face.size(); i++) {
			ave += verts[face[i] - 1];
		}
		ave /= (float)divisor;
		facePoints.push_back(ave);
	}
	// Then calculate edge points
	map<stuple<int, 2>, Vector3f> edgePs = edgePoints(edges, verts, facePoints, edgeFs);
	// And the average midpoints for each vertex
	vector<Vector3f> averageMidpts = aveMidEdges(verts, edges);
	// Next, calculate new vertices
	vector<Vector3f> newVerts;
	for (int i = 0; i < (int)verts.size(); i++) {
		int numFaces = (int)vertFs[i].size();
		Vector3f oldCoords = verts[i];
		Vector3f t1 = ((float) numFaces - 3) * oldCoords;
		Vector3f facePtAves = Vector3f(0.0f);
		for (int face : vertFs[i]) {
			facePtAves += facePoints[face];
		}
		facePtAves /= (float) numFaces;
		Vector3f aveEdges = averageMidpts[i];
		Vector3f newCoords = t1 + facePtAves + 2 * aveEdges;
		newCoords /= (float)numFaces;
		newVerts.push_back(newCoords);
	}
	// Finally, set up the faces
	vector<Vector3f> newPointsPlaceholder;
	vector<vector<int>> newFacesPlaceholder;
	vector<Vector3f>().swap(newPoints); // Clear and reallocate
	vector<vector<int>>().swap(newFaces);
	newPointsPlaceholder.insert(newPointsPlaceholder.end(), newVerts.begin(), newVerts.end());
	for (int f = 0; f < (int)faces.size(); f++) {
		vector<int> currentFace = faces[f];
		for (int v = 0; v < (int)currentFace.size(); v++) {
			int currentVert = currentFace[v];
			vector<int> newFace = { currentFace[v] };
			int nextVert;
			if (v + 1 == (int)currentFace.size()) {
				nextVert = currentFace[0];
			} else {
				nextVert = currentFace[v + 1];
			}
			int lastVert;
			if (v - 1 == -1) {
				lastVert = currentFace[(int)currentFace.size() - 1];
			} else {
				lastVert = currentFace[v - 1];
			}
			// Per point, new face has the vertices: new coordinate (already added), 
			// edge point from current point to next point in the original face, the face point,
			// and finally the edge point 
			// When each point is found, look it up and add it to the vertices if not already there
			stuple<int, 2> firstEdge = stuple<int, 2>(min(currentVert, nextVert), max(currentVert, nextVert));
			Vector3f firstEdgePoint = edgePs[firstEdge];
			bool isNew = true;
			int pos;
			for (int p = 0; p < (int)newPointsPlaceholder.size(); p++) {
				if (newPointsPlaceholder[p] == firstEdgePoint) {
					isNew = false;
					pos = p + 1;
				}
			}
			if (isNew) {
				newPointsPlaceholder.push_back(firstEdgePoint);
				newFace.push_back((int)newPointsPlaceholder.size());
			} else {
				newFace.push_back(pos);
			}
			Vector3f facePoint = facePoints[f]; 
			bool isNewF = true;
			int posF;
			for (int pf = 0; pf < (int)newPointsPlaceholder.size(); pf++) {
				if (newPointsPlaceholder[pf] == facePoint) {
					isNewF = false;
					posF = pf + 1;
				}
			}
			if (isNewF) {
				newPointsPlaceholder.push_back(facePoint);
				newFace.push_back((int)newPointsPlaceholder.size());
			} else {
				newFace.push_back(posF);
			}
			
			stuple<int, 2> lastEdge = stuple<int, 2>(min(currentVert, lastVert), max(currentVert, lastVert));
			Vector3f secondEdgePoint = edgePs[lastEdge];
			bool isNew2 = true;
			int pos2;
			for (int p2 = 0; p2 < (int)newPointsPlaceholder.size(); p2++) {
				if (newPointsPlaceholder[p2] == secondEdgePoint) {
					isNew2 = false;
					pos2 = p2 + 1;
				}
			}
			if (isNew2) {
				newPointsPlaceholder.push_back(secondEdgePoint);
				newFace.push_back((int)newPointsPlaceholder.size());
			} else {
				newFace.push_back(pos2);
			}
			newFacesPlaceholder.push_back(newFace);
		}
	}
	if (reps == 1) {
		newPoints.swap(newPointsPlaceholder);
		newFaces.swap(newFacesPlaceholder);
	} else {
		vector<Vector3f> newPtPlace2;
		vector<vector<int>> newFcPlace2;
		subsurfCC(newPointsPlaceholder, newFacesPlaceholder, newPtPlace2, newFcPlace2, reps - 1);
		newPoints.swap(newPtPlace2);
		newFaces.swap(newFcPlace2);
	}
}



// Main routine.
// Set up OpenGL, define the callbacks and start the main loop
int main(int argc, char** argv) {
	loadInput();

	GLFWwindow* window = createOpenGLWindow(600, 600, "Subsurf Test");


	// setup the keyboard event handler
	glfwSetKeyCallback(window, keyCallback);
	glfwSetMouseButtonCallback(window, mouseCallback);
	glfwSetCursorPosCallback(window, motionCallback);

	// glEnable() and glDisable() control parts of OpenGL's
	// fixed-function pipeline, such as rasterization, or
	// depth-buffering. What happens if you remove the next line?
	glClearColor(0, 0, 0, 1);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	// The program object controls the programmable parts
	// of OpenGL. All OpenGL programs define a vertex shader
	// and a fragment shader.
	program = compileProgram(c_vertexshader, c_fragmentshader);
	if (!program) {
		printf("Cannot compile program\n");
		return -1;
	}
	program_color = compileProgram(c_vertexshader_color, c_fragmentshader_color);
	if (!program_color) {
		printf("Cannot compile program\n");
		return -1;
	}

	camera.SetDimensions(600, 600);
	camera.SetPerspective(50);
	camera.SetDistance(10);
	camera.SetCenter(Vector3f(0, 0, 0));

	// Strip normals from vecf in preparation for subsurfing
	vector<vector<int>> strippedF;
	for (vector<int> face : vecf) {
		if ((int)face.size() != 4) {
			vector<int> newFace;
			for (int k = 0; k < (int)face.size(); k += 2) {
				newFace.push_back(face[k]);
			}
			strippedF.push_back(newFace);
		} else {
			strippedF.push_back(face);
		}
	}


	
	/*for (Vector3f vert : ssvecv) {
		cerr << "v " << vert[0] << " " << vert[1] << " " << vert[2] << endl;
	}
	for (vector<int> face : ssvecf) {
		cerr << "f";
		for (int p : face) {
			cerr << " " << p;
		}
		cerr << endl;
	}*/

    // Main Loop
    while (!glfwWindowShouldClose(window)) {
        // Clear the rendering window
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        setViewport(window);

		if (subsurfed > 0) {
			subsurfCC(vecv, strippedF, ssvecv, ssvecf, subsurfed);
		}
		if (!wireframe) {
			glUseProgram(program_color);
			camera.SetUniforms(program_color);
			updateLightUniforms(program_color);
			updateMaterialUniforms(program_color);
			updateRotationUniforms(program_color);
			glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		} else {
			glUseProgram(program);
			camera.SetUniforms(program);
			updateLightUniforms(program);
			updateMaterialUniforms(program);
			updateRotationUniforms(program);
			glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
			glLineWidth(1);
		}

        // Draw to back buffer
        drawScene();

        // Make back buffer visible
        glfwSwapBuffers(window);

        // Check if any input happened during the last frame
        glfwPollEvents();
    }

    // All OpenGL resource that are created with
    // glGen* or glCreate* must be freed.
    glDeleteProgram(program);

    glfwTerminate(); // destroy the window
    return 0;
}
