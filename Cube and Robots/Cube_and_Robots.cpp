#include "iostream"
#include "HW3.h"
using namespace std;
// OPENGL Variables
int th = 0;            //  Azimuth of view angle
int ph = 0;           //  Elevation of view angle
int axes = 1;         //  Display axes
int light = 1;
double asp = 1;     // aspect ratio
int fov = 45;         //  Field of view (for perspective)
double dim = 1.0;  // size of the world

//Simulator Variables
Cube cube[30];
std::vector<Gene> individual;

std::vector<std::vector<int>> population;
std::vector<std::vector<int>> selectedpop;

std::vector<double> popdistance;


double constantk = 3000;
double constantg = 9.81;
double constantr = 3000;
double constantf = 2.5;
double constants = 0.99;
double timestep = 0.001;
double T = 0;
int breath_flag = 0;
int mass_num = MASSNUM;
int spring_num = SPRINGNUM;
double potentialk = 0;
double potentialm = 0;
double potentialg = 0;
double potentialE = 0;
double kineticE = 0;
double totalE = 0;

std::ofstream bestgeneF("bestgene.txt");
std::ofstream dotplotF("dotplot.txt");
std::ofstream bestdistF("bestdist.txt");

GLfloat worldRotation[16] = { 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1 };


template <typename T>
vector<size_t> sort_indexes(const vector<T>& v) {

    // initialize original index locations
    vector<size_t> idx(v.size());
    iota(idx.begin(), idx.end(), 0);

    // sort indexes based on comparing values in v
    sort(idx.begin(), idx.end(),
        [&v](size_t i1, size_t i2) {return v[i1] < v[i2]; });

    return idx;
}

double normDistance(Mass m1,Mass m2) {
    double dx = m2.p[0] - m1.p[0];
    double dy = m2.p[1] - m1.p[1];
    double dz = m2.p[2] - m1.p[2];
    double dxy = std::sqrt(std::pow(dx, 2) + std::pow(dy, 2));
    double dxyz = std::sqrt(std::pow(dxy, 2) + std::pow(dz, 2));
    return dxyz;
}
void init_mass(double mass, double length, double X, double Y, double Z, int indexc) {
    Mass mass_base[MASSNUM];
    mass_base[0] = { mass,{X,Y,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[1] = { mass,{X + length,Y,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[2] = { mass,{X + 2 * length,Y,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[3] = { mass,{X,Y + length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[4] = { mass,{X + length,Y + length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[5] = { mass,{X + 2 * length,Y + length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[6] = { mass,{X,Y + 2 * length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[7] = { mass,{X + length,Y + 2 * length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[8] = { mass,{X + 2 * length,Y + 2 * length,Z},{0,0,0},{0,0,0},{0,0,0} };

    mass_base[9] = { mass,{X,Y,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[10] = { mass,{X + length,Y,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[11] = { mass,{X + 2 * length,Y,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[12] = { mass,{X,Y + length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[13] = { mass,{X + length,Y + length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[14] = { mass,{X + 2 * length,Y + length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[15] = { mass,{X,Y + 2 * length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[16] = { mass,{X + length,Y + 2 * length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[17] = { mass,{X + 2 * length,Y + 2 * length,Z + length},{0,0,0},{0,0,0},{0,0,0} };

    mass_base[18] = { mass,{X,Y,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[19] = { mass,{X + length,Y,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[20] = { mass,{X + 2 * length,Y,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[21] = { mass,{X,Y + length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[22] = { mass,{X + length,Y + length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[23] = { mass,{X + 2 * length,Y + length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[24] = { mass,{X,Y + 2 * length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[25] = { mass,{X + length,Y + 2 * length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[26] = { mass,{X + 2 * length,Y + 2 * length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };

    //srand((unsigned)time(NULL));
    vector<int> mass_element;
    for (int j=0; j < MASSNUM; j++) {
        double mass_flag = (double)(rand() / (double)RAND_MAX);
        if (mass_flag > 0.7) {
            mass_element.push_back(j);
        }
    }
    mass_num = int(mass_element.size());
    if (mass_num == 0) {
        mass_element.push_back(0);
        mass_element.push_back(2);
        mass_element.push_back(7);
        mass_element.push_back(9);
        mass_element.push_back(22);
    }
    cube[indexc].cube_mass.clear();
    cube[indexc].cube_spring.clear();
    mass_num = int(mass_element.size());
    for (int k = 0; k < mass_num; k++) {
        cube[indexc].cube_mass.push_back(mass_base[mass_element[k]]);
    }
    population.push_back(mass_element);
}

void init_spring(double constantk, int indexc) {
    Spring spring_temp;
    mass_num = int(cube[indexc].cube_mass.size());
    for (int i = 0; i < mass_num; i++) {
        for (int j = i + 1; j < mass_num; j++) {
            spring_temp = { 6000, normDistance(cube[indexc].cube_mass[i],cube[indexc].cube_mass[j]), normDistance(cube[indexc].cube_mass[i],cube[indexc].cube_mass[j]),i, j };
            cube[indexc].cube_spring.push_back(spring_temp);
        }
    }
    /*for (int i = 0; i < spring_num; i++) {
        int p1 = cube[indexc].cube_spring[i].m1;
        int p2 = cube[indexc].cube_spring[i].m2;
        cout << "p1: " << p1 << endl;
        cout << "p2: " << p2 << endl;
    }*/
    //cout << cube[indexc].cube_spring.size() << endl;
}

void init_cube(int cubenum) {
    
    double originalx = 0;
    double originaly = 0;
    double originalz = 0.2;
    
    for (int c = 0; c < cubenum; c++) {
        init_mass(0.1, 0.1, originalx, originaly, originalz, c);
        init_spring(constantk, c);
    }
}

void update_cube(int cubenum) {
    double originalx = 0;
    double originaly = 0;
    double originalz = 0.2;
    vector<int> mass_element;
    float mass = 0.1;
    float length = 0.1;
    float X = originalx;
    float Y = originaly;
    float Z = originalz;
    Mass mass_base[MASSNUM];
    mass_base[0] = { mass,{X,Y,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[1] = { mass,{X + length,Y,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[2] = { mass,{X + 2 * length,Y,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[3] = { mass,{X,Y + length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[4] = { mass,{X + length,Y + length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[5] = { mass,{X + 2 * length,Y + length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[6] = { mass,{X,Y + 2 * length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[7] = { mass,{X + length,Y + 2 * length,Z},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[8] = { mass,{X + 2 * length,Y + 2 * length,Z},{0,0,0},{0,0,0},{0,0,0} };

    mass_base[9] = { mass,{X,Y,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[10] = { mass,{X + length,Y,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[11] = { mass,{X + 2 * length,Y,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[12] = { mass,{X,Y + length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[13] = { mass,{X + length,Y + length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[14] = { mass,{X + 2 * length,Y + length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[15] = { mass,{X,Y + 2 * length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[16] = { mass,{X + length,Y + 2 * length,Z + length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[17] = { mass,{X + 2 * length,Y + 2 * length,Z + length},{0,0,0},{0,0,0},{0,0,0} };

    mass_base[18] = { mass,{X,Y,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[19] = { mass,{X + length,Y,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[20] = { mass,{X + 2 * length,Y,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[21] = { mass,{X,Y + length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[22] = { mass,{X + length,Y + length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[23] = { mass,{X + 2 * length,Y + length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[24] = { mass,{X,Y + 2 * length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[25] = { mass,{X + length,Y + 2 * length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    mass_base[26] = { mass,{X + 2 * length,Y + 2 * length,Z + 2 * length},{0,0,0},{0,0,0},{0,0,0} };
    
    for (int i = 0; i < cubenum; i++) {
        cube[i].cube_mass.clear();
        cube[i].cube_spring.clear();
        mass_element = population[i];
        mass_num = int(mass_element.size());
        for (int k = 0; k < mass_num; k++) {
            cube[i].cube_mass.push_back(mass_base[mass_element[k]]);
            bestgeneF << mass_element[k] << " ";
        }
        init_spring(constantk, i);
        bestgeneF << endl;
    }
    for (int c = cubenum; c < 2*cubenum; c++) {
        init_mass(0.1, 0.1, originalx, originaly, originalz, c);
        init_spring(constantk, c);
    }

}


void drawSomething()
{
    glColor3f(1, 0, 0);

    GLUquadric* quad;
    quad = gluNewQuadric();
    glPushMatrix();
    glMultMatrixf(worldRotation);
    glTranslated(0, 0, 0);
    gluSphere(quad, 1.0 / 10, 10, 10);
    glPopMatrix();

    glPushMatrix();
    glMultMatrixf(worldRotation);
    glBegin(GL_LINES);
    glVertex3f(0, 0, 0);
    glVertex3f(1, 1, 1);
    glEnd();
    glPopMatrix();

}

void forceCube(int cubenum) {
    T += timestep;
    //cout << T << endl;
    for (int c = 0; c < cubenum; c++) {
        mass_num = int(cube[c].cube_mass.size());
        spring_num = int(cube[c].cube_spring.size());
        potentialk = 0;
        potentialm = 0;
        potentialg = 0;
        potentialE = 0;
        kineticE = 0;
        totalE = 0;
        for (int j = 0; j < mass_num; j++) {
            cube[c].cube_mass[j].f[0] = 0;
            cube[c].cube_mass[j].f[1] = 0;
            cube[c].cube_mass[j].f[2] = -cube[c].cube_mass[j].m * constantg;
        }
        if (T >= 1.5){
            double geneb;
            double genec;
            for (int genum = 0; genum < spring_num; genum++) {
                int march_point = 1000;
                int ppp1 = cube[c].cube_spring[genum].m1;
                int ppp2 = cube[c].cube_spring[genum].m2;
                Mass massp1 = cube[c].cube_mass[ppp1];
                Mass massp2 = cube[c].cube_mass[ppp2];
                if (march_point > 0) {
                    if (massp1.p[2] < 0.001) {
                        march_point = ppp1;
                    }
                    if (massp2.p[2] < 0.001) {
                        march_point = ppp2;
                    }

                }
                if (ppp1 == march_point || ppp2 == march_point) {
                    geneb = 0.0225518;
                    genec = -4.17;
                    cube[c].cube_spring[genum].L0 = cube[c].cube_spring[genum].oriL + geneb / (2 * M_PI) * sin(10 * T + genec);
                }
                else {
                    geneb = 0.0116782;
                    genec = -3.86;
                    cube[c].cube_spring[genum].L0 = cube[c].cube_spring[genum].oriL + geneb / (2 * M_PI) * sin(10 * T + genec);
                }
                
            }
            /*for (int sn = 0; sn < SPRINGNUM; sn++) {
                    cube[c].cube_spring[sn].L0 = cube[c].cube_spring[sn].oriL + pow(-1, sn)* 0.2 * 0.005 * sin(10 * T);
                }*/
            
        }
        for (int i = 0; i < spring_num; i++) {
            int p1 = cube[c].cube_spring[i].m1;
            int p2 = cube[c].cube_spring[i].m2;
            //cout << "p1: " << p1 << endl;
            //cout << "p2: " << p2 << endl;
            Mass mass1 = cube[c].cube_mass[p1];
            Mass mass2 = cube[c].cube_mass[p2];
            double dxyz[3] = { mass2.p[0] - mass1.p[0],mass2.p[1] - mass1.p[1],mass2.p[2] - mass1.p[2] };
            double L = normDistance(mass1, mass2);
            double springforce = cube[c].cube_spring[i].k * (L - cube[c].cube_spring[i].L0);
        
            double forceDire[3] = { dxyz[0] / L,dxyz[1] / L,dxyz[2] / L };
            cube[c].cube_mass[p1].f[0] = cube[c].cube_mass[p1].f[0] + springforce * forceDire[0];
            cube[c].cube_mass[p1].f[1] = cube[c].cube_mass[p1].f[1] + springforce * forceDire[1];
            cube[c].cube_mass[p1].f[2] = cube[c].cube_mass[p1].f[2] + springforce * forceDire[2];
            cube[c].cube_mass[p2].f[0] = cube[c].cube_mass[p2].f[0] - springforce * forceDire[0];
            cube[c].cube_mass[p2].f[1] = cube[c].cube_mass[p2].f[1] - springforce * forceDire[1];
            cube[c].cube_mass[p2].f[2] = cube[c].cube_mass[p2].f[2] - springforce * forceDire[2];
            //potentialk = potentialk + cube[c].cube_spring[i].k * std::pow((L - cube[c].cube_spring[i].L0), 2) / 2;
        }
        
        for (int j = 0; j < mass_num; j++) {
            if (cube[c].cube_mass[j].p[2] < 0) {
                cube[c].cube_mass[j].f[2] = cube[c].cube_mass[j].f[2] + constantr * fabs(cube[c].cube_mass[j].p[2]);
                double hor_f = sqrt(pow(cube[c].cube_mass[j].f[0], 2) + pow(cube[c].cube_mass[j].f[1], 2));
                double ver_f = cube[c].cube_mass[j].f[2];
                if (hor_f < ver_f * constantf) {
                    cube[c].cube_mass[j].f[0] = 0;
                    cube[c].cube_mass[j].f[1] = 0;
                    cube[c].cube_mass[j].v[0] = 0;
                    cube[c].cube_mass[j].v[1] = 0;
                }
                else {
                    for (int h = 0; h < 2; h++) {
                        if (cube[c].cube_mass[j].f[h] > 0) {
                            cube[c].cube_mass[j].f[h] = cube[c].cube_mass[j].f[h] - ver_f * constantf * cube[c].cube_mass[j].f[h] / hor_f;
                            if (cube[c].cube_mass[j].f[h] < 0) {
                                cube[c].cube_mass[j].f[h] = 0;
                            }
                        }
                        else {
                            cube[c].cube_mass[j].f[h] = cube[c].cube_mass[j].f[h] + ver_f * constantf * fabs(cube[c].cube_mass[j].f[h]) / hor_f;
                            if (cube[c].cube_mass[j].f[h] > 0) {
                                cube[c].cube_mass[j].f[h] = 0;
                            }
                        }
                    }
                }
                //potentialg = potentialg + constantr * pow(cube[c].cube_mass[j].p[2], 2) / 2;
            }
            //potentialm = potentialm + cube[c].cube_mass[j].m * constantg * cube[c].cube_mass[j].p[2];
            //kineticE = kineticE + cube[c].cube_mass[j].m * (pow(cube[c].cube_mass[j].v[0], 2) + pow(cube[c].cube_mass[j].v[1], 2) + pow(cube[c].cube_mass[j].v[2], 2)) / 2;
            for (int a = 0; a < 3; a++) {
                cube[c].cube_mass[j].a[a] = cube[c].cube_mass[j].f[a] / cube[c].cube_mass[j].m;
                /*if (T < 0.0003) {
                    cout << "j: " << j << " a: " << a << cube[c].cube_mass[j].a[a] << endl;
                }*/
                //cout << "j: " << j << " a: " << a << cube[c].cube_mass[j].a[a] << endl;
            }
            for (int v = 0; v < 3; v++) {
                cube[c].cube_mass[j].v[v] = constants * (cube[c].cube_mass[j].v[v] + cube[c].cube_mass[j].a[v] * timestep);
                /*if (T < 0.0003) {
                    cout << "j: " << j << " v: " << v << cube[c].cube_mass[j].v[v] << endl;
                }*/
                //cout << "j: " << j << " v: " << v << cube[c].cube_mass[j].v[v] << endl;
            }
            for (int p = 0; p < 3; p++) {
                cube[c].cube_mass[j].p[p] = cube[c].cube_mass[j].p[p]+cube[c].cube_mass[j].v[p] * timestep;
                /*if (T < 0.0003) {
                    cout << "j: " << j << " p: " << p << cube[c].cube_mass[j].p[p] << endl;
                }*/
                //cout << "j: " << j << " p: " << p << cube[c].cube_mass[j].p[p] << endl;
            }

        }
        
        /*potentialE = potentialk + potentialg + potentialm;
        totalE = potentialE + kineticE;
        potentialF << T << " " << potentialE << endl;
        kineticF << T << " " << kineticE << endl;
        totalF << T << " " << totalE << endl;*/
        /*cout << T <<" " <<potentialE << endl;
        cout << T <<" " << kineticE << endl;
        cout << T << " " << totalE << endl;*/
    }
    
}
void drawground()
{
    glPushMatrix();
    glColor3f(0.96078, 0.96078, 0.86274);
    glBegin(GL_QUADS);
    glNormal3f(0, 1, 0);
    glTexCoord2f(0.0, 0.0);  glVertex3f(-2, 0, -2);
    glTexCoord2f(0.0, 1.0);  glVertex3f(-2, 0, +2);
    glTexCoord2f(1.0, 1.0);  glVertex3f(+2, 0, +2);
    glTexCoord2f(1.0, 0.0);  glVertex3f(+2, 0, -2);
    glEnd();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
    for (int i = 0; i < 10; i++) {
        for (int j = -9; j < 10; j++) {
            glColor3f(0, 1, 0);
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glBegin(GL_LINES);
            glLineWidth(10);
            glVertex3f(-0.5 * i / 10, 0.02, 0.5 * j / 10);
            glEnd();
            glPopMatrix();
        }
    }
}


void draw_single_spring(Spring cube_spring, int c) {
    if (c < 15) {
        int m1 = cube_spring.m1;
        int m2 = cube_spring.m2;
        glColor3f(0, 0, 1);
        glPushMatrix();
        glMultMatrixf(worldRotation);
        glBegin(GL_LINES);
        glVertex3f(cube[c].cube_mass[m1].p[0], cube[c].cube_mass[m1].p[1], cube[c].cube_mass[m1].p[2]);
        glVertex3f(cube[c].cube_mass[m2].p[0], cube[c].cube_mass[m2].p[1], cube[c].cube_mass[m2].p[2]);
        glEnd();
        glPopMatrix();
    }
    else {
        int m1 = cube_spring.m1;
        int m2 = cube_spring.m2;
        glColor3f(1, 0, 0);
        glPushMatrix();
        glMultMatrixf(worldRotation);
        glBegin(GL_LINES);
        glVertex3f(cube[c].cube_mass[m1].p[0], cube[c].cube_mass[m1].p[1], cube[c].cube_mass[m1].p[2]);
        glVertex3f(cube[c].cube_mass[m2].p[0], cube[c].cube_mass[m2].p[1], cube[c].cube_mass[m2].p[2]);
        glEnd();
        glPopMatrix();
    }
    
}


void drawCube( ) {
    if (T == 0) {
        std::cout << "population: " << population.size() << endl;
        best_gene_init();
        update_cube(population.size());
    }
    forceCube(30);
    if (T > 10) {
        cubeDistance(population.size());
        T = 0;
    }
    //drawground();
    /*for(int i=0;i<8;i++){
        std::cout << cube_mass[i].p[1] << std::endl;*/

    for (int c = 0; c < CUBENUM; c++) {
        if (c == 14) {
            glColor3f(1, 0, 0);
            mass_num = int(cube[c].cube_mass.size());
            GLUquadric* quad;
            for (int i = 0; i < mass_num; i++) {
                quad = gluNewQuadric();
                glPushMatrix();
                glMultMatrixf(worldRotation);
                glTranslatef(cube[c].cube_mass[i].p[0], cube[c].cube_mass[i].p[1], cube[c].cube_mass[i].p[2]);
                gluSphere(quad, 1.0 / 200, 10, 10);
                glPopMatrix();
            }
            spring_num = int(cube[c].cube_spring.size());
            for (int i = 0; i < spring_num; i++) {
                draw_single_spring(cube[c].cube_spring[i], c);
            }
        }
        /*if (c < 15) {
            glColor3f(1, 0, 0);
            mass_num = int(cube[c].cube_mass.size());
            GLUquadric* quad;
            for (int i = 0; i < mass_num; i++) {
                quad = gluNewQuadric();
                glPushMatrix();
                glMultMatrixf(worldRotation);
                glTranslatef(cube[c].cube_mass[i].p[0], cube[c].cube_mass[i].p[1], cube[c].cube_mass[i].p[2]);
                gluSphere(quad, 1.0 / 200, 10, 10);
                glPopMatrix();
            }
            spring_num = int(cube[c].cube_spring.size());
            for (int i = 0; i < spring_num; i++) {
                draw_single_spring(cube[c].cube_spring[i], c);
            }
        }
        else {
            glColor3f(0, 0, 1);
            mass_num = int(cube[c].cube_mass.size());
            GLUquadric* quad;
            for (int i = 0; i < mass_num; i++) {
                quad = gluNewQuadric();
                glPushMatrix();
                glMultMatrixf(worldRotation);
                glTranslatef(cube[c].cube_mass[i].p[0], cube[c].cube_mass[i].p[1], cube[c].cube_mass[i].p[2]);
                gluSphere(quad, 1.0 / 200, 10, 10);
                glPopMatrix();
            }
            spring_num = int(cube[c].cube_spring.size());
            for (int i = 0; i < spring_num; i++) {
                draw_single_spring(cube[c].cube_spring[i], c);
            }
        }*/
    }
    
}

//void init_population() {
//    srand((unsigned)time(NULL));
//    for (int i = 0; i < POPNUM; i++) {
//        std::vector<Gene> tempindividual;
//        for (int j = 0; j < CHROMENUM; j++) {
//            double geneb = (-2 * M_PI + (double)(rand() / (double)RAND_MAX) * 4 * M_PI)/500;
//            double genec = -2 * M_PI + (double)(rand() / (double)RAND_MAX) * 4 * M_PI;
//            Gene chrome{ geneb, genec };
//            tempindividual.push_back(chrome);
//        }
//        population.push_back(tempindividual);
//    }
//}

void cubeDistance(int popsize) {
    popdistance.clear();
    popdistance.shrink_to_fit();
    double indistance = 0;
    double maxdistance = 0;
    for (int i = 0; i < popsize; i++) {
        indistance = sqrt(pow(cube[i].cube_mass[0].p[0], 2) + pow(cube[i].cube_mass[0].p[1], 2));
        if (indistance > 100) {
            indistance = 0;
        }
        if (indistance < 0) {
            indistance = 0;
        }
        //cout << "indistance: "<<indistance << endl;
        maxdistance = std::fmax(maxdistance, indistance);
        //cout << "Maxdistance: " << maxdistance << endl;
        popdistance.push_back(indistance);
    }
    bestdistF << maxdistance << endl;
    for (int disnum = 0; disnum < popdistance.size(); disnum++) {
        dotplotF << popdistance[disnum]<< " ";
    }
    dotplotF << endl;

    //sort(popdistance.begin(), popdistance.end(), greater<float>());
    std::cout << "Maxdistance: "<< maxdistance << endl;
    std::vector<size_t> index = sort_indexes(popdistance);
    selectedpop.clear();
    selectedpop.shrink_to_fit();
    for (int i = popsize-1; i >= 0; i--) {
        //cout << "index: "<<index[i] << endl;
        selectedpop.push_back(population[index[i]]);
    }
    population.clear();
    population.shrink_to_fit();
    for (int i = 0; i < POPNUM; i++) {
        population.push_back(selectedpop[i]);
    }
}

//void cross_over() {
//    int popsize = population.size();
//    for (int i = 0; i < popsize / 2; i++) {
//        std::vector<double> parent1, parent2;
//        for (int j = 0; j < population[i].size(); j++) {
//            parent1.push_back(population[i][j].b);
//            parent1.push_back(population[i][j].c);
//            parent2.push_back(population[popsize-i-1][j].b);
//            parent2.push_back(population[popsize-i-1][j].c);
//        }
//
//        if (i == 0) {
//            for (int bestgenenum = 0; bestgenenum < parent1.size(); bestgenenum++) {
//                bestgeneF << parent1[bestgenenum] << " ";
//            }
//            bestgeneF << endl;
//        }
//
//        int point1 = floor((rand() / (double)RAND_MAX) * parent1.size() / 2) * 2;
//        int point2 = floor((rand() / (double)RAND_MAX) * parent1.size() / 2) * 2;
//        /*cout << "point1: " << point1 << endl;
//        cout << "point2: " << point2 << endl;*/
//        if (point2 < point1) {
//            int temp = point1;
//            point1 = point2;
//            point2 = temp;
//        }
//        if (point2 == point1) { point2 += 2; }
//        std::vector<double> child1, child2;
//        for (int crossnum = 0; crossnum < point1; crossnum++) {
//            child1.push_back(parent1[crossnum]);
//            child2.push_back(parent2[crossnum]);
//        }
//        for (int crossnum = point1; crossnum < point2; crossnum++) {
//            child1.push_back(parent2[crossnum]);
//            child2.push_back(parent1[crossnum]);
//        }
//        for (int crossnum = point2; crossnum < parent1.size(); crossnum++) {
//            child1.push_back(parent1[crossnum]);
//            child2.push_back(parent2[crossnum]);
//        }
//        child1 = mutation(child1);
//        child2 = mutation(child2);
//
//        Gene genetemp1, genetemp2;
//        std::vector<Gene> childpop1, childpop2;
//        for (int genenum = 0; genenum < child1.size(); genenum += 2) {
//            genetemp1.b = child1[genenum];
//            genetemp1.c = child1[genenum + 1];
//            childpop1.push_back(genetemp1);
//            genetemp2.b = child2[genenum];
//            genetemp2.c = child2[genenum + 1];
//            childpop2.push_back(genetemp2);
//        }
//        population.push_back(childpop1);
//        population.push_back(childpop2);
//    }
//}

//std::vector<double> mutation(std::vector<double> child) {
//    for (int i = 0; i < child.size(); i++) {
//        double mutationpro = rand() / (double)RAND_MAX;
//        //cout << "mutationpro: " << mutationpro << endl;
//        if(mutationpro>0.8){
//            if (i % 2 == 0) {
//                child[i] = (-2 * M_PI + (double)(rand() / (double)RAND_MAX) * 4 * M_PI) / 500;
//            }
//            else {
//                child[i] = -2 * M_PI + (double)(rand() / (double)RAND_MAX) * 4 * M_PI;
//            }
//        }
//    }
//    return child;
//}


void best_gene_init() {
    population.clear();
    vector<int> mass_element;
    mass_element = { 10, 13, 14, 15, 18, 19, 20, 22, 23 };
    population.push_back(mass_element);
    mass_element = { 0, 1, 9, 15, 16, 17, 19, 24, 25, 26 };
    population.push_back(mass_element);
    mass_element = { 0, 1, 9, 13, 17, 18, 19, 22, 24, 25 };
    population.push_back(mass_element);
    mass_element = { 0, 1, 2, 4, 5, 6, 11, 12, 15, 16, 19, 21, 24, 26 };
    population.push_back(mass_element);
    mass_element = { 0, 12, 16, 18, 19, 21, 24, 25, 26 };
    population.push_back(mass_element);
    mass_element = { 7, 8, 13, 17, 21, 23 };
    population.push_back(mass_element);
    mass_element = { 0, 3, 4, 7, 9, 10, 13, 14, 15, 16 };
    population.push_back(mass_element);
    mass_element = { 3, 6, 7, 8, 9, 14, 15, 18, 19, 22, 23 };
    population.push_back(mass_element);
    mass_element = { 1, 2, 3, 13, 14, 15, 20, 21, 22 };
    population.push_back(mass_element);
    mass_element = { 2, 5, 7, 11, 16, 18, 24, 25, 26 };
    population.push_back(mass_element);
    mass_element = { 5, 6, 7, 10, 11, 12, 15, 19, 22, 23, 25, 26 };
    population.push_back(mass_element);
    mass_element = { 6, 8, 17, 19, 23, 25 };
    population.push_back(mass_element);
    mass_element = { 0, 3, 6, 11, 15, 16, 17, 24 };
    population.push_back(mass_element);
    mass_element = { 2, 6, 15, 16, 19, 22 };
    population.push_back(mass_element);
    mass_element = { 2, 3, 4, 11, 12, 13, 23, 24, 25 };
    population.push_back(mass_element);
}

void Print(const char* format, ...)
{
    char    buf[LEN];
    char* ch = buf;
    va_list args;
    //  Turn the parameters into a character string
    va_start(args, format);
    vsnprintf(buf, LEN, format, args);
    va_end(args);
    //  Display the characters one at a time at the current raster position
    while (*ch)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *ch++);
}

/*
 *  OpenGL (GLUT) calls this routine to display the scene
 */
void display()
{
    const double len = 1.5;  //  Length of axes
    //  Erase the window and the depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //  Enable Z-buffering in OpenGL
    glEnable(GL_DEPTH_TEST);
    //  Undo previous transformations
    glLoadIdentity();
    //  Eye position
    double Ex = -2 * dim * Sin(th) * Cos(ph);
    double Ey = +2 * dim * Sin(ph);
    double Ez = +2 * dim * Cos(th) * Cos(ph);
    gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, Cos(ph), 0);


    //drawSomething();
    drawCube();

    //  Draw axes
    glColor3f(1, 1, 1);
    if (axes)
    {
        glBegin(GL_LINES);
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(len, 0.0, 0.0);
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(0.0, len, 0.0);
        glVertex3d(0.0, 0.0, 0.0);
        glVertex3d(0.0, 0.0, len);
        glEnd();
        //  Label axes
        glRasterPos3d(len, 0.0, 0.0);
        Print("X");
        glRasterPos3d(0.0, len, 0.0);
        Print("Z");
        glRasterPos3d(0.0, 0.0, len);
        Print("Y");
    }
    //  Render the scene
    glFlush();
    //  Make the rendered scene visible
    glutSwapBuffers();
}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
void special(int key, int x, int y)
{
    //  Right arrow key - increase angle by 5 degrees
    if (key == GLUT_KEY_RIGHT)
        th += 5;
    //  Left arrow key - decrease angle by 5 degrees
    else if (key == GLUT_KEY_LEFT)
        th -= 5;
    //  Up arrow key - increase elevation by 5 degrees
    else if (key == GLUT_KEY_UP)
    {
        if (ph + 5 < 90)
        {
            ph += 5;
        }
    }
    //  Down arrow key - decrease elevation by 5 degrees
    else if (key == GLUT_KEY_DOWN)
    {
        if (ph - 5 > 0)
        {
            ph -= 5;
        }
    }
    //  Keep angles to +/-360 degrees
    th %= 360;
    ph %= 360;
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  Set projection
 */
void Project(double fov, double asp, double dim)
{
    //  Tell OpenGL we want to manipulate the projection matrix
    glMatrixMode(GL_PROJECTION);
    //  Undo previous transformations
    glLoadIdentity();
    //  Perspective transformation
    if (fov)
        gluPerspective(fov, asp, dim / 16, 16 * dim);
    //  Orthogonal transformation
    else
        glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);
    //  Switch to manipulating the model matrix
    glMatrixMode(GL_MODELVIEW);
    //  Undo previous transformations
    glLoadIdentity();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
void key(unsigned char ch, int x, int y)
{
    //  Exit on ESC
    if (ch == 27)
        exit(0);
    //  Reset view angle
    else if (ch == '0')
        th = ph = 0;
    //  Toggle axes
    else if (ch == 'a' || ch == 'A')
        axes = 1 - axes;
    //  Change field of view angle
    else if (ch == '-' && ch > 1)
        fov++;
    else if (ch == '=' && ch < 179)
        fov--;
    //  PageUp key - increase dim
    else if (ch == GLUT_KEY_PAGE_DOWN) {
        dim += 0.1;
    }
    //  PageDown key - decrease dim
    else if (ch == GLUT_KEY_PAGE_UP && dim > 1) {
        dim -= 0.1;
    }
    //  Keep angles to +/-360 degrees
    th %= 360;
    ph %= 360;
    //  Reproject
    Project(fov, asp, dim);
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  GLUT calls this routine when the window is resized
 */
void reshape(int width, int height)
{
    //  Ratio of the width to the height of the window
    asp = (height > 0) ? (double)width / height : 1;
    //  Set the viewport to the entire window
    glViewport(0, 0, width, height);
    //  Set projection
    Project(fov, asp, dim);
}

/*
 *  GLUT calls this toutine when there is nothing else to do
 */
void idle()
{
    glutPostRedisplay();
}

int main(int argc, char* argv[])
{
    init_cube(15);
    //init_population();
    //best_gene_init();
    // Initialize GLUT and process user parameters
    glutInit(&argc, argv);
    // double buffered, true color 600*600
    glutInitWindowSize(1000, 800);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    // create the window
    glutCreateWindow("xm2225");
    //  Tell GLUT to call "idle" when there is nothing else to do
    glutIdleFunc(idle);
    //  Tell GLUT to call "display" when the scene should be drawn
    glutDisplayFunc(display);
    //  Tell GLUT to call "reshape" when the window is resized
    glutReshapeFunc(reshape);
    //  Tell GLUT to call "special" when an arrow key is pressed
    glutSpecialFunc(special);
    //  Tell GLUT to call "key" when a key is pressed
    glutKeyboardFunc(key);
    //  Pass control to GLUT so it can interact with the user
    glutMainLoop();
    return 0;
    
    //init_cube(15);
    //while (1) {
    //    cout << "population: "<<population.size() << endl;
    //    update_cube(population.size());
    //    /*for (int i = 0; i < population.size(); i++) {
    //        for (int j = 0; j < population[i].size(); j++) {
    //            cout << population[i][j];
    //        }
    //        cout << " " << endl;
    //    }*/
    //    while (T < 10) {
    //        forceCube(population.size());
    //        T += timestep;
    //    }
    //    //cout << "i'm on the ground" << endl;
    //    cubeDistance(population.size());
    //    T = 0;
    //}
};


//#include "HW3.h"
//GLfloat worldRotation[16] = { 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1 };
//
//struct MASS
//{
//    double m;       // mass
//    double p[3];    // 3D position
//    double v[3];    // 3D velocity
//    double a[3];    // 3D acceleration
//};
//
//struct SPRING
//{
//    double k;       // spring constant
//    double L_0;     // rest length
//    int m1;         // first mass connected
//    int m2;         // second mass connected
//};
//
//struct GENE {
//    double k;
//    double b;
//    double c;
//};
////int num_foot = 2 + rand() % 6;
//
////ofstream outFile1("distance3.txt");
////ofstream outFile2("bestgene3.txt");
////ofstream outFile3("dot3.txt");
//
////struct Cube {
////    struct MASS cubemass[8];
////    struct SPRING cubespring[28];
////};
//
//float L(MASS mass1, MASS mass2) {
//    double length = sqrt(pow((mass1.p[0] - mass2.p[0]), 2) + pow((mass1.p[1] - mass2.p[1]), 2) + pow((mass1.p[2] - mass2.p[2]), 2));
//
//    return length;
//}
//
//
////vector<MASS> joint = jointmass(mass, length, 0, 0, 0.01);
////
////
////vector<SPRING> spring = cubespring(length, k);
//
//
//GLuint tex;
//GLUquadric* sphere;
//void make_tex(void)
//{
//    unsigned char data[256][256][3];
//    for (int y = 0; y < 255; y++) {
//        for (int x = 0; x < 255; x++) {
//            unsigned char* p = data[y][x];
//            p[0] = p[1] = p[2] = (x ^ y) & 8 ? 255 : 0;
//        }
//    }
//    glGenTextures(1, &tex);
//    glBindTexture(GL_TEXTURE_2D, tex);
//    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, (const GLvoid*)data);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//}
//
//void init(void)
//{
//    glEnable(GL_DEPTH_TEST);
//    make_tex();
//    sphere = gluNewQuadric();
//    glEnable(GL_TEXTURE_2D);
//}
//
//
//
//class Robot {
//private:
//    vector<GENE> walkgene;
//
//public:
//    vector<MASS> joint;
//    vector<SPRING> spring;
//    int mass_num;
//    int num_sp;
//    int num_foot;
//    vector<double> original;
//    double initl[3] = { 0,0,0 };
//    Robot(double x, double y, double z, int num, vector<GENE> newGene)
//    {
//        mass_num = 7 + num;
//        num_foot = num;
//        num_sp = 21+ 6 * num_foot;
//        initl[0] = x;
//        initl[1] = y;
//        initl[2] = z;
//        walkgene = newGene;
//        joint = jointmass(x, y, z, num);
//        cubespring(walkgene);
//
//    }
//
//    vector<MASS> jointmass(double x, double y, double z, int num)
//    {
//        vector<MASS> joint(7 + num);
//        double subtle = length / 2;
//        joint[0] = { mass, {x + 1.5 * length,y + length * sqrt(3) / 2, z + 1.5 * length}, {0,0,0}, {0,0,0} };
//        joint[1] = { mass, {x + length   ,y + sqrt(3) * length,z + 0.5 * length}, {0,0,0}, {0,0,0} };
//        joint[2] = { mass, {x + 2 * length  ,y + sqrt(3) * length ,z + 0.5 * length}, {0,0,0}, {0,0,0} };
//        joint[3] = { mass, {x + length / 2,y + length / 2 * sqrt(3),z + 0.5 * length}, {0,0,0}, {0,0,0} };
//        joint[4] = { mass, {x + length / 2 + 2 * length ,y + length / 2 * sqrt(3),z + 0.5 * length}, {0,0,0}, {0,0,0} };
//        joint[5] = { mass, {x + length,y,z + 0.5 * length}, {0,0,0}, {0,0,0} };
//        joint[6] = { mass, {x + 2 * length ,y ,z + 0.5 * length}, {0,0,0}, {0,0,0} };
//        for (int i = 7; i < 7 + num; i++) {
//            double a = -length + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 * length)));
//            double b = -length + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 * length)));
//            joint[i] = { mass, {joint[0].p[0] + a ,joint[0].p[1] + b, z}, {0,0,0}, {0,0,0} };
//        }
//        //joint[7] = { mass, {x + 1.5 * length ,y + length * 3 * sqrt(3) / 2 - subtle,z}, {0,0,0}, {0,0,0} };
//        //joint[8] = { mass, {x + subtle ,y + sqrt(3) * length - subtle ,z}, {0,0,0}, {0,0,0} };
//        //joint[9] = { mass, {x + 3 * length - subtle,y + sqrt(3) * length - subtle,z }, {0,0,0}, {0,0,0} };
//        //joint[10] = { mass, {x + subtle ,y + subtle ,z }, {0,0,0}, {0,0,0} };
//        //joint[11] = { mass, {x + 3 * length - subtle, y + subtle ,z}, {0,0,0}, {0,0,0} };
//        //joint[12] = { mass, {x + 1.5 * length,y - length * sqrt(3) / 2 + subtle, z}, {0,0,0}, {0,0,0} };
//        return joint;
//    }
//    void cubespring(vector<GENE> walkgene)
//    {
//
//        int pointer = 0;
//        for (int i = 0; i < 6; i++) {
//            for (int j = i + 1; j < 7; j++) {
//                original.push_back(L(joint[i], joint[j]));
//                spring.push_back({ walkgene[pointer].k,L(joint[i],joint[j]),i,j });
////                cout << spring[i * j].k << endl;
//                pointer++;
//
//
//            }
//        }
//        for (int i = 7; i < 7 + num_foot; i++) {
//            for (int j = 0; j < 6; j++) {
//                original.push_back(L(joint[i], joint[j]));
//                spring.push_back({ walkgene[pointer].k,L(joint[i],joint[j]),i,j });
//                pointer++;
//            }
//        }
//
//        //     spring[0] = { walkgene[0].k,L(joint[0],joint[3]),0,3 };
//        //     spring[1] = { walkgene[1].k,L(joint[0],joint[4]),0,4 };
//        //     spring[2] = { walkgene[2].k,L(joint[0],joint[5]),0,5 };
//        //
//        //     spring[3] = { walkgene[3].k,L(joint[1],joint[3]),1,3 };
//        //     spring[4] = { walkgene[4].k,L(joint[1],joint[5]),1,5 };
//        //     spring[5] = { walkgene[5].k,L(joint[1],joint[4]),1,4 };
//        //
//        //     spring[6] = { walkgene[6].k,L(joint[2],joint[4]),2,4 };
//        //     spring[7] = { walkgene[7].k,L(joint[2],joint[5]),2,5 };
//        //     spring[8] = { walkgene[8].k,L(joint[2],joint[3]),2,3 };
//        //
//        //     spring[9] = { walkgene[9].k,L(joint[4],joint[5]),4,5 };
//        //     spring[10] = { walkgene[10].k,L(joint[4],joint[3]),4,3 };
//        //     spring[11] = { walkgene[11].k,L(joint[5],joint[3]),3,5 };
//        //
//        //     spring[12] = { walkgene[12].k,L(joint[3],joint[6]),3,6 };
//        //     spring[13] = { walkgene[13].k,L(joint[4],joint[6]),4,6 };
//        //     spring[14] = { walkgene[14].k,L(joint[5],joint[6]),5,6 };
//        //
//        //     spring[15] = { walkgene[15].k,L(joint[0],joint[6]),0,6 };
//        //     spring[16] = { walkgene[16].k,L(joint[1],joint[6]),1,6 };
//        //     spring[17] = { walkgene[17].k,L(joint[2],joint[6]),2,6 };
//        //
//        //     spring[18] = { walkgene[18].k,L(joint[3],joint[7]),3,7 };
//        //     spring[19] = { walkgene[19].k,L(joint[4],joint[7]),4,7 };
//        //     spring[20] = { walkgene[20].k,L(joint[5],joint[7]),5,7 };
//        //        spring[21] = { walkgene[21].k,L(joint[6],joint[7]),6,7 };
//        //
//        //    spring[22] = { walkgene[22].k,L(joint[3],joint[8]),3,8 };
//        //    spring[23] = { walkgene[23].k,L(joint[4],joint[8]),4,8 };
//        //    spring[24] = { walkgene[24].k,L(joint[5],joint[8]),5,8 };
//        //    spring[25] = { walkgene[25].k,L(joint[6],joint[8]),6,8 };
//
//    }
//
//    void drawcube()
//    {
//
//        glPushMatrix();
//        glMultMatrixf(worldRotation);
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(R,G,B);
//        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
//        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
//        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(R,G,B);
//        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
//        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
//        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(R,G,B);
//        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
//        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
//        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(R,G,B);
//        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
//        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
//        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(R,G,B);
//        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
//        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
//        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(R,G,B);
//        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
//        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
//        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
//        glEnd();
//        //foot
////        glBegin(GL_TRIANGLES);
////        glColor3f(0.2, 0.1, 0.3);
////        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
////        glVertex3f(GLfloat(joint[7].p[0]), GLfloat(joint[7].p[1]), GLfloat(joint[7].p[2]));
////        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
////        glEnd();
////
////        glBegin(GL_TRIANGLES);
////        glColor3f(0.2, 0.1, 0.3);
////        glVertex3f(GLfloat(joint[8].p[0]), GLfloat(joint[8].p[1]), GLfloat(joint[8].p[2]));
////        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
////        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
////        glEnd();
////
////        glBegin(GL_TRIANGLES);
////        glColor3f(0.2, 0.1, 0.3);
////        glVertex3f(GLfloat(joint[9].p[0]), GLfloat(joint[9].p[1]), GLfloat(joint[9].p[2]));
////        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
////        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
////        glEnd();
////
////        glBegin(GL_TRIANGLES);
////        glColor3f(0.2, 0.1, 0.3);
////        glVertex3f(GLfloat(joint[11].p[0]), GLfloat(joint[11].p[1]), GLfloat(joint[11].p[2]));
////        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
////        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
////        glEnd();
////
////        glBegin(GL_TRIANGLES);
////        glColor3f(0.2, 0.1, 0.3);
////        glVertex3f(GLfloat(joint[12].p[0]), GLfloat(joint[12].p[1]), GLfloat(joint[12].p[2]));
////        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
////        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
////        glEnd();
////
////        glBegin(GL_TRIANGLES);
////        glColor3f(0.2, 0.1, 0.3);
////        glVertex3f(GLfloat(joint[10].p[0]), GLfloat(joint[10].p[1]), GLfloat(joint[10].p[2]));
////        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
////        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
////        glEnd();
////
//        glPopMatrix();
//
//
//        //
//        GLUquadric* quad[mass_num];
//        for (int i = 1; i < mass_num; i++) {
//
//            glColor3f(R1,G1, B1);
//            quad[i] = gluNewQuadric();
//            glPushMatrix();
//            glMultMatrixf(worldRotation);
//            glTranslated(joint[i].p[0], joint[i].p[1], joint[i].p[2]);
//
//            gluSphere(quad[i], 1.0 / 250, 10, 10);
//            //            if (i == 6) { gluSphere(quad[i], 1.0 / 100, 10, 10); }
//            glPopMatrix();
//
//        }
//
//        for (int i = 0; i < num_sp; i++) {
//            int a = spring[i].m1;
//            int b = spring[i].m2;
//
//            glColor3f(0.5, 0.5, 0.5);
//            glPushMatrix();
//            glMultMatrixf(worldRotation);
//            glBegin(GL_LINES);
//            glLineWidth(20);
//            glVertex3f(joint[a].p[0], joint[a].p[1], joint[a].p[2]);
//            glVertex3f(joint[b].p[0], joint[b].p[1], joint[b].p[2]);
//            glEnd();
//            glPopMatrix();
//
//        }
//    }
//
//    void drawground()
//    {
//        glPushMatrix();
//        glColor3f(0.96078, 0.96078, 0.86274);
//        glBegin(GL_QUADS);
//        glNormal3f(0, 1, 0);
//        glTexCoord2f(0.0, 0.0);  glVertex3f(-1.5, +0.0, -1.5);
//        glTexCoord2f(0.0, 1);  glVertex3f(+1.5, +0.0, -1.5);
//        glTexCoord2f(1, 1);  glVertex3f(+1.5, +0.0, +1.5);
//        glTexCoord2f(1, 0.0);  glVertex3f(-1.5, +0.0, +1.5);
//        glEnd();
//        glPopMatrix();
//        glDisable(GL_TEXTURE_2D);
//        for (int i = 0; i < 10; i++) {
//            for (int j = -9; j < 10; j++) {
//                glColor3f(0, 1, 0);
//                glPushMatrix();
//                glMultMatrixf(worldRotation);
//                glBegin(GL_LINES);
//                glLineWidth(10);
//                glVertex3f(-0.5 * i / 10, 0.02, 0.5 * j / 10);
//                //        glV/Users/chiqu/class2019fall/EA/assignment3/HW3/HW3.cppertex3f(0.5*i/10, 0.5*j/10, 0.01);
//                glEnd();
//                glPopMatrix();
//            }
//        }
//    }
//
//    void simulate() {
//        vector<vector<double>> cForces(num_foot + 7, vector<double>(3));
//        for (int i = 0; i < mass_num; i++) {
//            cForces[i][0] = 0.0;
//            cForces[i][1] = 0.0;
//            cForces[i][2] = -joint[i].m * gravity;
//        }
//
//        for (int i = 0; i < num_sp; i++) {
//            MASS mass1 = joint[spring[i].m1];
//            MASS mass2 = joint[spring[i].m2];
//            //if (T < 1) {
////            if (i == 18 or i == 37 or i == 17 or i == 27 or i == 29 or i == 39 or i == 46 or i == 48 or i == 54 or i == 56 or i == 62 or i == 64) {
////
////            spring[i].L_0 = original[i] * (1 + walkgene[i].b * sin(500 * T + walkgene[i].c));
////            if (i > num_sp / 2) {
////                spring[i].L_0 = original[i] * (1 + walkgene[i].b * cos(500 * T + walkgene[i].c));
////            }}
//            if (i > 21 && i < 21+ 6 * num_foot) {
//                spring[i].L_0 = original[i] * (1 + walkgene[i].b * sin(500 * T + walkgene[i].c));
//            }
//
//            double pd[3] = { mass2.p[0] - mass1.p[0],mass2.p[1] - mass1.p[1],mass2.p[2] - mass1.p[2] };
//            double new_L = L(mass1, mass2);
//            double L_0 = spring[i].L_0;
//            double force = spring[i].k * fabs(new_L - L_0);
//            //springenergy += k * pow((new_L - L_0), 2) / 2;
//            //cout <<i<<"---new_L---" <<new_L << endl;
//            double norm_pd[3] = { pd[0] / new_L, pd[1] / new_L, pd[2] / new_L };
//            //cout << i << "---force---" << force << endl;
//            //compression
//            if (new_L < spring[i].L_0) {
//                cForces[spring[i].m1][0] -= norm_pd[0] * force;
//                cForces[spring[i].m1][1] -= norm_pd[1] * force;
//                cForces[spring[i].m1][2] -= norm_pd[2] * force;
//                cForces[spring[i].m2][0] += norm_pd[0] * force;
//                cForces[spring[i].m2][1] += norm_pd[1] * force;
//                cForces[spring[i].m2][2] += norm_pd[2] * force;
//            }
//
//            //tension
//            else {
//                cForces[spring[i].m1][0] += norm_pd[0] * force;
//                cForces[spring[i].m1][1] += norm_pd[1] * force;
//                cForces[spring[i].m1][2] += norm_pd[2] * force;
//                cForces[spring[i].m2][0] -= norm_pd[0] * force;
//                cForces[spring[i].m2][1] -= norm_pd[1] * force;
//                cForces[spring[i].m2][2] -= norm_pd[2] * force;
//            }
//        }
//        //cout << T << "tan " << springenergy << endl;
//        for (int i = 0; i < mass_num; i++) {
//            if (joint[i].p[2] <= 0) {
//                cForces[i][2] -= Nground * joint[i].p[2];
//                //cForces[i][2] -= joint[i].m*0.981;
//                //groundenergy += Nground * pow(joint[i].p[2], 2) / 2;
//                double Fh = sqrt(pow(cForces[i][0], 2) + pow(cForces[i][1], 2));
//                double Fv = cForces[i][2];
//                if (Fh < Fv * frictionCoefficient) {
//                    cForces[i][0] = 0;
//                    cForces[i][1] = 0;
//                    joint[i].v[0] = 0;
//                    joint[i].v[1] = 0;
//                }
//                else {
//                    double Fh_new = Fh - Fv * frictionCoefficient;
//                    cForces[i][0] = cForces[i][0] * Fh_new / Fh;
//                    cForces[i][1] = cForces[i][1] * Fh_new / Fh;
//                }
//            }
//            //cout << T << " " << groundenergy << endl;
//            for (int j = 0; j < 3; j++) {
//                joint[i].a[j] = cForces[i][j] / joint[i].m;
//                joint[i].v[j] += joint[i].a[j] * timeStep;
//                joint[i].v[j] *= dampening;
//                joint[i].p[j] += joint[i].v[j] * timeStep;
//            }
//            //cout << cube[i].p[j] << endl;
//        //gravityenergy
//            //gravityenergy += joint[i].m * gravity * joint[i].p[2];
//            //kinetic
//            double norm_v = sqrt(pow(joint[i].v[0], 2) + pow(joint[i].v[1], 2) + pow(joint[i].v[2], 2));
//            //kineticenergy += joint[i].m * pow(norm_v, 2) / 2;
//
//
//            //            cout << i << joint[i].p[2] << endl;
//            //            cout << cForces[1][2] << endl;
//        }
//        //cout << "hhhhh"<< endl;
//
//        //cout << "wtfh" << endl;
//
//        drawground();
//
//        T = T + timeStep;
//
//    }
//
//};
//
//vector<vector<GENE>> populationGene;
//vector<Robot> robot;
//vector<double> totald2;
//vector<vector<GENE>> newgenepop;
//vector<vector<GENE>> nextpop;
//vector<GENE> bestGene;
//int populationsize = size;
//
//
//
//vector<double> bestdistance;
//
////double maxdistance;
//vector<vector<GENE>> bestgene;
//
//
//vector<GENE> crossover(vector<GENE> p1, vector<GENE> p2) {
//    if (p1.size() < p2.size()) {
//        int crossposition = 7 + rand() % static_cast<int>(p1.size());
//
//        vector<GENE> offspring = p1;
//        //cout << "crossover" << endl;
//        for (int i = 0; i < crossposition; i++) {
//            offspring[i] = p2[i];
//        }
//
//        return offspring;
//    }
//
//    else {
//        int crossposition = 7 + rand() % static_cast<int>(p2.size());
//        vector<GENE> offspring = p2;
//        //cout << "crossover" << endl;
//        for (int i = 0; i < crossposition; i++) {
//            offspring[i] = p1[i];
//        }
//
//        return offspring;
//    }
//}
//
//vector<vector<GENE>> genereateGene(int num_foot) {
//    vector<vector<GENE>> populationGene1;
//    for (int i = 0; i < populationsize; i++) {
//        vector<GENE> temp;
//        GENE temp_2;
////        vector<GENE> temp_4 = {{9000,0,0},{900,0, 0},{9000,0, 0},{900,0, 0},{9000, 0, 0},{9000, 0, 0},{9000,0, 0},{9000, 0, 0},{9000, 0, 0},{9000, 0, 0},{9000, 0, 0},{9000, 0, 0},
////        {9000, 0, 0},
////        {9000, 0, 0},
////        {9000, 0, 0},
////        {9000, 0, 0},
////        {9000, 0, 0},
////        {9000, 0, 0},
////        {9000, 0, 0},
////        {9000,0, 0},
////        {9000, 0, 0},
////        {9000, 0, 0},
////        {1149.03, 0.89871, 1.59599},
////        {3387.09, 0.0564407, 1.23423},
////        {2002.62, 0.699834, -4.80743},
////        {2914.94, 0.553911, 1.1057},
////        {817.105, 0.321366, -3.82508},
////        {4385, 0.813266, 0.71422},
////        {1138.03, 0.402719, -0.0315316},
////        {896.933, 0.0647815, 3.55796},
////        {1613.19, 0.143684, 5.02053},
////        {1201.31, 0.177645, 2.19434},
////        {682.43, 0.771142, 1.13884},
////        {1140.98, 0.416475, 2.54433},
////        {1074.05, 0.766175,-5.02824},
////        {2116.14, 0.297064, 3.23736},
////        {2055.34, 0.0556062, 0.919644},
////        {794.846, 0.637075, -2.20964},
////        {2764.1, 0.367395, 3.83164},
////        {1421.09, 0.473704, 0.641933},
////        {2274.68, 0.883705, -0.953366},
////        {1192.37, 0.17468, 4.34046},
////        {1291.49, 0.273336, 5.86253},
////        {1497.38, 0.0790328, -2.46469},
////        {2897.75, 0.861622, -2.79981},
////        {932.343, 0.309218, -5.87749},
////        {1368.84, 0.606336, 2.27619},
////        {963.681, 0.530644, 0.467611},
////        {870.933, 0.522875, 5.70392},
////        {4467.87, 0.359379, -5.15406},
////            {993.408, 0.668242, -4.40994}};
//        double k_1 = 9000;
//        double b_1 = 0;
//        double c_1 = 0;
//        int num_sp = 21 + 6 * num_foot;
//        for (int j = 0; j < num_sp; j++) {
////            double k_1 = 100 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (1000 - 100)));
////            //                if (j == 18 or j == 37 or j == 17 or j == 27 or j == 29 or j == 39 or j == 46 or j == 48 or j == 54 or j == 56 or j == 62 or j == 64) {
////            if (j < num_sp && j > 21) {
////                k_1 = 600 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (5000 - 600)));
////                b_1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 1);
////                c_1 = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));
////            }
////            else {
////                k_1 = 9000;
////                b_1 = 0;
////                c_1 = 0;
//            k_1 = 9000;
//            b_1 = 0;
//            c_1 = 0;
//
//
//            temp_2 = { k_1, b_1, c_1 };
//
//            temp.push_back(temp_2);}
//            //                tempVec.push_back(k1);
//            //                tempVec.push_back(b1);
//            //                tempVec.push_back(c1);
//
//
//        populationGene1.push_back(temp);
//
//    }
//
//    return populationGene1;
//}
//
//vector<Robot> sortrobot(vector<Robot> robot, vector<double>totald){
//
////     vector<int> Index;
////         for (int i = 0; i < populationsize; i++) {
////             cout<<"totald "<<totald[i]<<endl;
////             double temp_min = totald[i];
////             int pointer = i;
////             for (int j = i; j < populationsize; j++) {
////                 if (totald[j] > temp_min) {
////                     temp_min = totald[j];
////                     pointer = j;
////
////                 }
////
////             }
////             sort(totald.begin(),totald.end());
////             double a = *max_element (totald.begin(), totald.end());
////             k = distance(totald,a);
////             cout<<"pointer"<<pointer<<endl;
////             Index.push_back(pointer);
////         }
//        vector<size_t> idx(totald.size());
//        iota(idx.begin(), idx.end(), 0);
//
//        // sort indexes based on comparing values in v
//        sort(idx.begin(), idx.end(),
//             [totald](size_t i1, size_t i2) {return totald[i1] > totald[i2];});
//
//    vector<Robot> tempRobot;
//    for(int i = 0; i<populationsize;i++){
//
//
//
//        Robot tempp = robot[idx[i]];
//        tempRobot.push_back(tempp);
//    }
//
//    return tempRobot;
//
//}
//
//
//
//
//
//vector<vector<GENE>> selection(vector<vector<GENE>> populationGene, vector<double> totald) {
//    vector<vector<GENE>> newgenepop;
//    vector<size_t> idx(totald.size());
//    iota(idx.begin(), idx.end(), 0);
//
//    // sort indexes based on comparing values in v
//    sort(idx.begin(), idx.end(),
//         [totald](size_t i1, size_t i2) {return totald[i1] > totald[i2];});
//
//   vector<size_t> Index = idx;
////
////        for (int i = 0; i < populationsize; i++) {
////            double temp_min = totald[0];
////            int pointer = 0;
////            for (int j = i; j < populationsize; j++) {
////                if (totald[j] > temp_min) {
////                    temp_min = totald[j];
////                    pointer = j;
////                }
////
////            }
////            Index.push_back(pointer);
////        }
//
//    newgenepop.clear();
//    newgenepop.shrink_to_fit();
//    for (int i = 0; i < Index.size() / 2; i++) {
//        newgenepop.push_back(populationGene[Index[i]]);
//    }
//    bestgene.push_back(populationGene[Index[0]]);
//    bestdistance.push_back(totald[Index[0]]);
//    cout << "best distance" << totald[Index[0]] << endl;
//    evaluation++;
//    //
////    outFile1<<evaluation<<" "<< totald[Index[0]]<<endl;
//
//    cout << "evaluation: " << evaluation << endl;
//            //cout<<"worst"<<totald[index[-1]]<<endl;
////    for (int i=0;i<newgenepop[0].size(); i++){
////        outFile2<<i <<" "<< newgenepop[0][i].k << " " <<newgenepop[0][i].b<< " " <<newgenepop[0][i].c << " ";
////    }
////
//    return newgenepop;
//}
//vector<GENE> genereateNewGene(int num_foot) {
//
//        vector<GENE> temp;
//    GENE temp_2;
//
////        double k_1 = 9000;
////        double b_1 = 0;
////        double c_1 = 0;
////        int num_sp = 21+ 6 * num_foot;
////
////        for (int j = 0; j < num_sp; j++) {
////            double k_1 = 100 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (1000 - 100)));
////            //                if (j == 18 or j == 37 or j == 17 or j == 27 or j == 29 or j == 39 or j == 46 or j == 48 or j == 54 or j == 56 or j == 62 or j == 64) {
////            if (j < num_sp && j > 21) {
////                k_1 = 600 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (5000 - 600)));
////                b_1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 1);
////                c_1 = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));
////            }
////            else {
////                k_1 = 9000;
////                b_1 = 0;
////                c_1 = 0;
////            }
//
////            temp_2 = { k_1, b_1, c_1 };
//            temp.push_back(temp_2);
//
//
//
//    return temp;
//}
//Robot getnewrobot(vector<GENE> gene, int num_foot){
//
////double x = -2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 4));
//double y = -1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2));
////    cout<<"getnewgenenenene"<<gene[9].k<<endl;
//Robot robot2(0, y, 0.1, num_foot, gene);
//
//
//
//return robot2;
//}
//
//vector<GENE> mutation(vector<GENE> p1) {
//    vector<int> list;
//    for (int i = 21; i < p1.size(); i++) {
//        list.push_back(i);
//    }
//    int i = 21 + rand() % static_cast<int>(p1.size() - 21);
//    int i1 = list[i];
//    vector<GENE> offspring = p1;
//    GENE temp;
//    temp.k = 600 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (1000 - 100)));
//    temp.b = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 1);
//    temp.c = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));
//
//    offspring[i] = temp;
//    list.clear();
//    return offspring;
//}
//vector<vector<GENE>> nextgeneration(vector<Robot> robots,vector<vector<GENE>> newgenepop) {
//    vector<vector<GENE>> nextpop = newgenepop;
//    vector<GENE> p3;
//    for (int a = 0; a<int(newgenepop.size())-1; a++) {
//        int i = rand() % static_cast<int>(newgenepop.size());
//        int j = rand() % static_cast<int>(newgenepop.size());
//
//        vector<GENE> p1 = nextpop[i];
//        vector<GENE> p2 = newgenepop[j];
//        if (p1.size() < p2.size()) {
//            p1 = crossover(p1, p2);
//            nextpop.push_back(p1);
//        }
//        else {
//            p2 = crossover(p1, p2);
//            nextpop.push_back(p2);
//        }
//
//    }
//    int num = 2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 10));
//
//
//
//    vector<GENE> temp3 = genereateNewGene(num);
//    nextpop.push_back(temp3);
//
//
//
//    robots.pop_back();
//    //nextpop[-1] = temp3;
////     cout <<"nextpop[-1].size(2222)" << nextpop[-1][5].k <<  endl;
//
//    Robot a = getnewrobot(nextpop[9], num);
//    robots.push_back(a);
//
////    robots[-1] = getnewrobot(nextpop[-1], num);
//    double r = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
//    if (r < 0.5) {
//        vector<GENE> temp_off;
//        int ra_ro = 1 + rand() % static_cast<int>(nextpop.size() - 1);
//        temp_off = mutation(nextpop[ra_ro]);
//        nextpop[ra_ro] = temp_off;
//        cout << "mutationing" << endl;
//
//
//    }
//            for(int i =0; i<nextpop[0].size(); i++){
////                outFile2<<evaluation<<" "<<i <<" "<< nextpop[0][i].k<<" "<<nextpop[0][i].b<<" "<<nextpop[0][i].c<<endl;
//            }
//    return nextpop;
//}
//
////vector<Robot> getrobots(Robot a){
////    robots.pop_back();
////    robots.push_back(a);
////    return robots;
////}
////
//
//
//
//
//vector<Robot> getinitalrobot(vector<vector<GENE>> nextpop, int num_foot) {
//    vector<Robot> robots;
//    for (int i = 0; i < populationsize; i++) {
//
//        //double x = -2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 4));
//        double y = -1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2));
//        //robots.push_back(ROBOT(0.0, 3.0*(i-populationSize/2), 0.0, populationGene[i]));
////            int num = 2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 8));
//
//        Robot robot1(0, y, 0.1, num_foot, nextpop[i]);
//        robots.push_back(robot1);
//
//    }
//    return robots;
//}
//
//    double fittness(Robot r, int num_m) {
//       // double fit = 0;
//        double startpoint;
//        double sumdis = 0;
//        for (int i = 0; i < num_m-1; i++) {
////            fit += (r.joint[i].p[0] - r.initl[0]);
//            sumdis += (r.joint[i].p[0] - r.initl[0]);
//
//        }
//        startpoint = sumdis / num_m;
//        r.initl[0] = startpoint;
////        outFile3<< evaluation << " " << sumdis/(num_m-1) << endl;
//                //cout<<"startpoint"<<startpoint<<endl;
//        return sumdis / (num_m-1);
//    }
//
//
//    vector<double> totaldistance(vector<Robot> robots) {
//        vector<double> totald;
//        for (int i = 0; i < populationsize; i++) {
//            double f = fittness(robots[i], robots[i].num_foot + 7);
//
//            totald.push_back(f);
//
//        }
//
//        return totald;
//    }
//
//
//
//
//
//
//
//    void Print(const char* format, ...)
//    {
//        char    buf[LEN];
//        char* ch = buf;
//        va_list args;
//        //  Turn the parameters into a character string
//        va_start(args, format);
//        vsnprintf(buf, LEN, format, args);
//        va_end(args);
//        //  Display the characters one at a time at the current raster position
//        while (*ch)
//            glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *ch++);
//    }
//
//    /*
//     *  OpenGL (GLUT) calls this routine to display the scene
//     */
//
//
//    void display()
//    {
//        //drawground();
//        double len = 0.2;  //  Length of axes
//        //  Erase the window and the de5th buffer
//        glClearColor(0.5372549, 0.6549019, 0.760784, 1.0);
//        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//        //  Enable Z-buffering in OpenGL
//        glEnable(GL_DEPTH_TEST);
//        //  Undo previous transformations
//        glLoadIdentity();
//        //  Eye position
//        double Ex = -1 * dim * Sin(th) * Cos(ph);
//        double Ey = +1 * dim * Sin(ph);
//        double Ez = +1 * dim * Cos(th) * Cos(ph);
//        gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, Cos(ph), 0);
//
//        //Simulate();
//
////        double y = -1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2));
////
////        Robot robot1(0, y, 0.1, 5, populationGene);
////
////        if (T >0) {
//
//            for (int j = 0; j < populationsize; j++) {
////                cout<<robot[j].num_foot<<endl;
//                if(j==0){
//                robot[j].simulate();
//
//                robot[j].drawcube();
//                }
////            }
//
//
//        }
////        else {
////            totald2 = totaldistance(robot);
////            newgenepop = selection(populationGene, totald2);
////            nextpop = nextgeneration(robot,newgenepop);
////            robot = sortrobot(robot, totald2);
////            robot.pop_back();
////            robot.push_back(getnewrobot(nextpop[9], int (nextpop[9].size()-21)/6));
////            populationGene = nextpop;
////            //cout<<"in loop evolve"<<endl;
////            T = 0;
////        }
//
//        //Draw axes
//        //glColor3f(1, 0, 0);
//    //    if (axes)
//    //    {
//    //        glBegin(GL_LINES);
//    //        glLineWidth(2);
//    //        glVertex3d(0.0, 0.0, 0.0);
//    //        glVertex3d(len, 0.0, 0.0);
//    //        glVertex3d(0.0, 0.0, 0.0);
//    //        glVertex3d(0.0, len, 0.0);
//    //        glVertex3d(0.0, 0.0, 0.0);
//    //        glVertex3d(0.0, 0.0, len);
//    //        glEnd();
//    //        //  Label axes
//    //        glRasterPos3d(len, 0.0, 0.0);
//    //        Print("X");
//    //        glRasterPos3d(0.0, len, 0.0);
//    //        Print("Y");
//    //        glRasterPos3d(0.0, 0.0, len);
//    //        Print("Z");
//    //    }
//        //  Render the scene
//        glFlush();
//        //  Make the rendered scene visible
//        glutSwapBuffers();
//    }
//
//
//    void special(int key, int x, int y)
//    {
//        //  Right arrow key - increase angle by 5 degrees
//        if (key == GLUT_KEY_RIGHT)
//            th += 5;
//        //  Left arrow key - decrease angle by 5 degrees
//        else if (key == GLUT_KEY_LEFT)
//            th -= 5;
//        //  Up arrow key - increase elevation by 5 degrees
//        else if (key == GLUT_KEY_UP)
//        {
//            if (ph + 5 < 90)
//            {
//                ph += 5;
//            }
//        }
//        //  Down arrow key - decrease elevation by 5 degrees
//        else if (key == GLUT_KEY_DOWN)
//        {
//            if (ph - 5 > 0)
//            {
//                ph -= 5;
//            }
//        }
//        //  Keep angles to +/-360 degrees
//        th %= 360;
//        ph %= 360;
//        //  Tell GLUT it is necessary to redisplay the scene
//        glutPostRedisplay();
//    }
//
//
//    void Project(double fov, double asp, double dim)
//    {
//        //  Tell OpenGL we want to manipulate the projection matrix
//        glMatrixMode(GL_PROJECTION);
//        //  Undo previous transformations
//        glLoadIdentity();
//        //  Perspective transformation
//        if (fov)
//            gluPerspective(fov, asp, dim / 16, 16 * dim);
//        //  Orthogonal transformation
//        else
//            glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);
//        //  Switch to manipulating the model matrix
//        glMatrixMode(GL_MODELVIEW);
//        //  Undo previous transformations
//        glLoadIdentity();
//    }
//
//
//    void key(unsigned char ch, int x, int y)
//    {
//        //  Exit on ESC
//        if (ch == 27)
//            exit(0);
//        //  Reset view angle
//        else if (ch == '0')
//            th = ph = 0;
//        //  Toggle axes
//        else if (ch == 'a' || ch == 'A')
//            axes = 1 - axes;
//        //  Change field of view angle
//        else if (ch == '-' && ch > 1)
//            fov++;
//        else if (ch == '=' && ch < 179)
//            fov--;
//        //  PageUp key - increase dim
//        else if (ch == GLUT_KEY_PAGE_DOWN) {
//            dim += 0.1;
//        }
//        //  PageDown key - decrease dim
//        else if (ch == GLUT_KEY_PAGE_UP && dim > 1) {
//            dim -= 0.1;
//        }
//        //  Keep angles to +/-360 degrees
//        th %= 360;
//        ph %= 360;
//        //  Reproject
//        Project(fov, asp, dim);
//        //  Tell GLUT it is necessary to redisplay the scene
//        glutPostRedisplay();
//    }
//
//
//    void reshape(int width, int height)
//    {
//        //  Ratio of the width to the height of the window
//        asp = (height > 0) ? (double)width / height : 1;
//        //  Set the viewport to the entire window
//        glViewport(0, 0, width, height);
//        //  Set projection
//        Project(fov, asp, dim);
//    }
//
//
//    void idle()
//    {
//        glutPostRedisplay();
//    }
//int foot = 5;
//
//
//    int main(int argc, char* argv[])
//    {
//        populationGene = genereateGene(foot);
//
//        robot = getinitalrobot(populationGene,5);
////        Robot robot1 = getnewrobot(populationGene, 5);
//
//        // Initialize GLUT and process user parameters
//        glutInit(&argc, argv);
//        // double buffered, true color 600*600
//        glutInitWindowSize(1000, 800);
//        glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
//        // create the window
//        glutCreateWindow("evolve");
//        //  Tell GLUT to call "idle" when there is nothing else to do
//        glutIdleFunc(idle);
//        //  Tell GLUT to call "display" when the scene should be drawn
//        glutDisplayFunc(display);
//        //  Tell GLUT to call "reshape" when the window is resized
//        glutReshapeFunc(reshape);
//        //  Tell GLUT to call "special" when an arrow key is pressed
//        glutSpecialFunc(special);
//        //  Tell GLUT to call "key" when a key is pressed
//        glutKeyboardFunc(key);
//        init();
//
//        //  Pass control to GLUT so it can interact with the user
//        glutMainLoop();
//
//        return 0;
//
//
//
//    };
//hw5 ends here

/*
#include "HW3.h"
//#include "Header.h"
GLfloat worldRotation[16] = { 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1 };

struct MASS
{
    double m;       // mass
    double p[3];    // 3D position
    double v[3];    // 3D velocity
    double a[3];    // 3D acceleration
};

struct SPRING
{
    double k;       // spring constant
    double L_0;     // rest length
    int m1;         // first mass connected
    int m2;         // second mass connected
};

struct GENE {
    double k;
    double b;
    double c;
};
//int num_foot = 2 + rand() % 6;

//ofstream outFile1("distance.txt");
//ofstream outFile2("bestgene.txt");
//ofstream outFile3("dot.txt");

//struct Cube {
//    struct MASS cubemass[8];
//    struct SPRING cubespring[28];
//};

float L(MASS mass1, MASS mass2) {
    double length = sqrt(pow((mass1.p[0] - mass2.p[0]), 2) + pow((mass1.p[1] - mass2.p[1]), 2) + pow((mass1.p[2] - mass2.p[2]), 2));

    return length;
}


//vector<MASS> joint = jointmass(mass, length, 0, 0, 0.01);
//
//
//vector<SPRING> spring = cubespring(length, k);


GLuint tex;
GLUquadric* sphere;
void make_tex(void)
{
    unsigned char data[256][256][3];
    for (int y = 0; y < 255; y++) {
        for (int x = 0; x < 255; x++) {
            unsigned char* p = data[y][x];
            p[0] = p[1] = p[2] = (x ^ y) & 8 ? 255 : 0;
        }
    }
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, (const GLvoid*)data);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
}

void init(void)
{
    glEnable(GL_DEPTH_TEST);
    make_tex();
    sphere = gluNewQuadric();
    glEnable(GL_TEXTURE_2D);
}



class Robot {
private:
    vector<GENE> walkgene;

public:
    vector<MASS> joint;
    vector<SPRING> spring;
    int mass_num;
    int num_sp;
    int num_foot;
    vector<double> original;
    double initl[3] = { 0,0,0 };
    Robot(double x, double y, double z, int num, vector<GENE> newGene)
    {
        mass_num = 7+num;
        num_foot = num;
        num_sp = 15+6*num_foot;
        initl[0] = x;
        initl[1] = y;
        initl[2] = z;
        walkgene = newGene;
        joint = jointmass(x, y, z, num);
        cubespring(walkgene);
        
    }

    vector<MASS> jointmass(double x, double y, double z, int num)
    {
     vector<MASS> joint(7+num);
     double subtle = length / 2;
     joint[0] = { mass, {x + 1.5 * length,y + length * sqrt(3) / 2, z + 1.5 * length}, {0,0,0}, {0,0,0} };
     joint[1] = { mass, {x + length   ,y + sqrt(3) * length,z + 0.5 * length}, {0,0,0}, {0,0,0} };
     joint[2] = { mass, {x + 2 * length  ,y + sqrt(3) * length ,z + 0.5 * length}, {0,0,0}, {0,0,0} };
     joint[3] = { mass, {x + length / 2,y + length / 2 * sqrt(3),z + 0.5 * length}, {0,0,0}, {0,0,0} };
     joint[4] = { mass, {x + length / 2 + 2 * length ,y + length / 2 * sqrt(3),z + 0.5 * length}, {0,0,0}, {0,0,0} };
     joint[5] = { mass, {x + length,y,z + 0.5 * length}, {0,0,0}, {0,0,0} };
     joint[6] = { mass, {x + 2 * length / 2,y + length / (2 * sqrt(3)),z + 0.5 * length}, {0,0,0}, {0,0,0} };
     for (int i = 7; i < 7+num; i++) {
      double a = -length + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 * length)));
      double b = -length + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (2 * length)));
      joint[i] = { mass, {joint[0].p[0] + a ,joint[0].p[1] + b, z}, {0,0,0}, {0,0,0} };
     }
     //joint[7] = { mass, {x + 1.5 * length ,y + length * 3 * sqrt(3) / 2 - subtle,z}, {0,0,0}, {0,0,0} };
     //joint[8] = { mass, {x + subtle ,y + sqrt(3) * length - subtle ,z}, {0,0,0}, {0,0,0} };
     //joint[9] = { mass, {x + 3 * length - subtle,y + sqrt(3) * length - subtle,z }, {0,0,0}, {0,0,0} };
     //joint[10] = { mass, {x + subtle ,y + subtle ,z }, {0,0,0}, {0,0,0} };
     //joint[11] = { mass, {x + 3 * length - subtle, y + subtle ,z}, {0,0,0}, {0,0,0} };
     //joint[12] = { mass, {x + 1.5 * length,y - length * sqrt(3) / 2 + subtle, z}, {0,0,0}, {0,0,0} };
     return joint;
    }
    void cubespring(vector<GENE> walkgene)
    {

        int pointer = 0;
        for (int i = 0; i < 6; i++) {
            for (int j = i + 1; j < 7; j++) {
                original.push_back(L(joint[i],joint[j]));
                spring.push_back( { walkgene[pointer].k,L(joint[i],joint[j]),i,j });
                cout<<spring[i*j].k<<endl;
                pointer++;


            }
        }
        for (int i = 7; i < 7+num_foot; i++){
            for (int j = 0; j < 6 ; j++){
                original.push_back(L(joint[i],joint[j]));
                spring.push_back( { walkgene[pointer].k,L(joint[i],joint[j]),i,j });
                pointer++;
            }
        }
        
        //     spring[0] = { walkgene[0].k,L(joint[0],joint[3]),0,3 };
        //     spring[1] = { walkgene[1].k,L(joint[0],joint[4]),0,4 };
        //     spring[2] = { walkgene[2].k,L(joint[0],joint[5]),0,5 };
        //
        //     spring[3] = { walkgene[3].k,L(joint[1],joint[3]),1,3 };
        //     spring[4] = { walkgene[4].k,L(joint[1],joint[5]),1,5 };
        //     spring[5] = { walkgene[5].k,L(joint[1],joint[4]),1,4 };
        //
        //     spring[6] = { walkgene[6].k,L(joint[2],joint[4]),2,4 };
        //     spring[7] = { walkgene[7].k,L(joint[2],joint[5]),2,5 };
        //     spring[8] = { walkgene[8].k,L(joint[2],joint[3]),2,3 };
        //
        //     spring[9] = { walkgene[9].k,L(joint[4],joint[5]),4,5 };
        //     spring[10] = { walkgene[10].k,L(joint[4],joint[3]),4,3 };
        //     spring[11] = { walkgene[11].k,L(joint[5],joint[3]),3,5 };
        //
        //     spring[12] = { walkgene[12].k,L(joint[3],joint[6]),3,6 };
        //     spring[13] = { walkgene[13].k,L(joint[4],joint[6]),4,6 };
        //     spring[14] = { walkgene[14].k,L(joint[5],joint[6]),5,6 };
        //
        //     spring[15] = { walkgene[15].k,L(joint[0],joint[6]),0,6 };
        //     spring[16] = { walkgene[16].k,L(joint[1],joint[6]),1,6 };
        //     spring[17] = { walkgene[17].k,L(joint[2],joint[6]),2,6 };
        //
        //     spring[18] = { walkgene[18].k,L(joint[3],joint[7]),3,7 };
        //     spring[19] = { walkgene[19].k,L(joint[4],joint[7]),4,7 };
        //     spring[20] = { walkgene[20].k,L(joint[5],joint[7]),5,7 };
        //        spring[21] = { walkgene[21].k,L(joint[6],joint[7]),6,7 };
        //
        //    spring[22] = { walkgene[22].k,L(joint[3],joint[8]),3,8 };
        //    spring[23] = { walkgene[23].k,L(joint[4],joint[8]),4,8 };
        //    spring[24] = { walkgene[24].k,L(joint[5],joint[8]),5,8 };
        //    spring[25] = { walkgene[25].k,L(joint[6],joint[8]),6,8 };
        cout<< spring[0].k <<endl;
        cout<<"creatspring"<<endl;
        
    }
    
    void drawcube()
    {
   
        glPushMatrix();
        glMultMatrixf(worldRotation);

        glBegin(GL_TRIANGLES);
        glColor3f(240.0/255, 136.0/255, 56.0/255);
        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
        glEnd();

        glBegin(GL_TRIANGLES);
        glColor3f(240.0/255, 136.0/255, 56.0/255);
        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
        glEnd();

        glBegin(GL_TRIANGLES);
        glColor3f(240.0/255, 136.0/255, 56.0/255);
        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
        glEnd();

        glBegin(GL_TRIANGLES);
        glColor3f(240.0/255, 136.0/255, 56.0/255);
        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
        glEnd();

        glBegin(GL_TRIANGLES);
        glColor3f(240.0/255, 136.0/255, 56.0/255);
        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
        glEnd();

        glBegin(GL_TRIANGLES);
        glColor3f(240.0/255, 136.0/255, 56.0/255);
        glVertex3f(GLfloat(joint[0].p[0]), GLfloat(joint[0].p[1]), GLfloat(joint[0].p[2]));
        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
        glEnd();
        //foot
//        glBegin(GL_TRIANGLES);
//        glColor3f(0.2, 0.1, 0.3);
//        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
//        glVertex3f(GLfloat(joint[7].p[0]), GLfloat(joint[7].p[1]), GLfloat(joint[7].p[2]));
//        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(0.2, 0.1, 0.3);
//        glVertex3f(GLfloat(joint[8].p[0]), GLfloat(joint[8].p[1]), GLfloat(joint[8].p[2]));
//        glVertex3f(GLfloat(joint[1].p[0]), GLfloat(joint[1].p[1]), GLfloat(joint[1].p[2]));
//        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(0.2, 0.1, 0.3);
//        glVertex3f(GLfloat(joint[9].p[0]), GLfloat(joint[9].p[1]), GLfloat(joint[9].p[2]));
//        glVertex3f(GLfloat(joint[2].p[0]), GLfloat(joint[2].p[1]), GLfloat(joint[2].p[2]));
//        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(0.2, 0.1, 0.3);
//        glVertex3f(GLfloat(joint[11].p[0]), GLfloat(joint[11].p[1]), GLfloat(joint[11].p[2]));
//        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
//        glVertex3f(GLfloat(joint[4].p[0]), GLfloat(joint[4].p[1]), GLfloat(joint[4].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(0.2, 0.1, 0.3);
//        glVertex3f(GLfloat(joint[12].p[0]), GLfloat(joint[12].p[1]), GLfloat(joint[12].p[2]));
//        glVertex3f(GLfloat(joint[6].p[0]), GLfloat(joint[6].p[1]), GLfloat(joint[6].p[2]));
//        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
//        glEnd();
//
//        glBegin(GL_TRIANGLES);
//        glColor3f(0.2, 0.1, 0.3);
//        glVertex3f(GLfloat(joint[10].p[0]), GLfloat(joint[10].p[1]), GLfloat(joint[10].p[2]));
//        glVertex3f(GLfloat(joint[3].p[0]), GLfloat(joint[3].p[1]), GLfloat(joint[3].p[2]));
//        glVertex3f(GLfloat(joint[5].p[0]), GLfloat(joint[5].p[1]), GLfloat(joint[5].p[2]));
//        glEnd();
//
        glPopMatrix();

        
//
        GLUquadric* quad[mass_num];
        for (int i = 0; i < mass_num; i++) {

            glColor3f(0.2, 0.4, 0.4);
            if (i == 6) { glColor3f(1, 0, 0); }
            quad[i] = gluNewQuadric();
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glTranslated(joint[i].p[0], joint[i].p[1], joint[i].p[2]);

            gluSphere(quad[i], 1.0 / 200, 10, 10);
//            if (i == 6) { gluSphere(quad[i], 1.0 / 100, 10, 10); }
            glPopMatrix();

        }

        for (int i = 0; i < num_sp; i++) {
            int a = spring[i].m1;
            int b = spring[i].m2;

            glColor3f(0.5, 0.5, 0.5);
            glPushMatrix();
            glMultMatrixf(worldRotation);
            glBegin(GL_LINES);
            glLineWidth(20);
            glVertex3f(joint[a].p[0], joint[a].p[1], joint[a].p[2]);
            glVertex3f(joint[b].p[0], joint[b].p[1], joint[b].p[2]);
            glEnd();
            glPopMatrix();

        }
    }

    void drawground()
    {
        glPushMatrix();
        glColor3f(0.96078, 0.96078, 0.86274);
        glBegin(GL_QUADS);
        glNormal3f(0, 1, 0);
        glTexCoord2f(0.0, 0.0);  glVertex3f(-1.5, +0.0, -1.5);
        glTexCoord2f(0.0, 1);  glVertex3f(+1.5, +0.0, -1.5);
        glTexCoord2f(1, 1);  glVertex3f(+1.5, +0.0, +1.5);
        glTexCoord2f(1, 0.0);  glVertex3f(-1.5, +0.0, +1.5);
        glEnd();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);
        for (int i = 0; i < 10; i++) {
            for (int j = -9; j < 10; j++) {
                glColor3f(0, 1, 0);
                glPushMatrix();
                glMultMatrixf(worldRotation);
                glBegin(GL_LINES);
                glLineWidth(10);
                glVertex3f(-0.5 * i / 10, 0.02, 0.5 * j / 10);
                //        glV/Users/chiqu/class2019fall/EA/assignment3/HW3/HW3.cppertex3f(0.5*i/10, 0.5*j/10, 0.01);
                glEnd();
                glPopMatrix();
            }
        }
    }

    void simulate() {
        vector<vector<double>> cForces(num_foot+7, vector<double>(3));
        for (int i = 0; i < mass_num; i++) {
            cForces[i][0] = 0.0;
            cForces[i][1] = 0.0;
            if (T == 1 or T == 10) {
                cForces[i][0] = 10.0;
            }

            cForces[i][2] = -joint[i].m * gravity;
        }

        for (int i = 0; i < num_sp; i++) {
            MASS mass1 = joint[spring[i].m1];
            MASS mass2 = joint[spring[i].m2];
            //if (T < 1) {
//            if (i == 18 or i == 37 or i == 17 or i == 27 or i == 29 or i == 39 or i == 46 or i == 48 or i == 54 or i == 56 or i == 62 or i == 64) {
//
//            spring[i].L_0 = original[i] * (1 + walkgene[i].b * sin(500 * T + walkgene[i].c));
//            if (i > num_sp / 2) {
//                spring[i].L_0 = original[i] * (1 + walkgene[i].b * cos(500 * T + walkgene[i].c));
//            }}
            if (i > 15 && i< 15+6*num_foot){
                spring[i].L_0 = original[i] * (1 + walkgene[i].b * sin(500 * T + walkgene[i].c));
            }
                
            double pd[3] = { mass2.p[0] - mass1.p[0],mass2.p[1] - mass1.p[1],mass2.p[2] - mass1.p[2] };
            double new_L = L(mass1, mass2);
            double L_0 = spring[i].L_0;
            double force = spring[i].k * fabs(new_L - L_0);
            springenergy += k * pow((new_L - L_0), 2) / 2;
            //cout <<i<<"---new_L---" <<new_L << endl;
            double norm_pd[3] = { pd[0] / new_L, pd[1] / new_L, pd[2] / new_L };
            //cout << i << "---force---" << force << endl;
            //compression
            if (new_L < spring[i].L_0) {
                cForces[spring[i].m1][0] -= norm_pd[0] * force;
                cForces[spring[i].m1][1] -= norm_pd[1] * force;
                cForces[spring[i].m1][2] -= norm_pd[2] * force;
                cForces[spring[i].m2][0] += norm_pd[0] * force;
                cForces[spring[i].m2][1] += norm_pd[1] * force;
                cForces[spring[i].m2][2] += norm_pd[2] * force;
            }

            //tension
            else {
                cForces[spring[i].m1][0] += norm_pd[0] * force;
                cForces[spring[i].m1][1] += norm_pd[1] * force;
                cForces[spring[i].m1][2] += norm_pd[2] * force;
                cForces[spring[i].m2][0] -= norm_pd[0] * force;
                cForces[spring[i].m2][1] -= norm_pd[1] * force;
                cForces[spring[i].m2][2] -= norm_pd[2] * force;
            }
        }
        //cout << T << "tan " << springenergy << endl;
        for (int i = 0; i < mass_num; i++) {
            if (joint[i].p[2] <= 0) {
                cForces[i][2] -= Nground * joint[i].p[2];
                //cForces[i][2] -= joint[i].m*0.981;
                groundenergy += Nground * pow(joint[i].p[2], 2) / 2;
                double Fh = sqrt(pow(cForces[i][0], 2) + pow(cForces[i][1], 2));
                double Fv = cForces[i][2];
                if (Fh < Fv * frictionCoefficient) {
                    cForces[i][0] = 0;
                    cForces[i][1] = 0;
                    joint[i].v[0] = 0;
                    joint[i].v[1] = 0;
                }
                else {
                    double Fh_new = Fh - Fv * frictionCoefficient;
                    cForces[i][0] = cForces[i][0] * Fh_new / Fh;
                    cForces[i][1] = cForces[i][1] * Fh_new / Fh;
                }
            }
            //cout << T << " " << groundenergy << endl;
            for (int j = 0; j < 3; j++) {
                joint[i].a[j] = cForces[i][j] / joint[i].m;
                joint[i].v[j] += joint[i].a[j] * timeStep;
                joint[i].v[j] *= dampening;
                joint[i].p[j] += joint[i].v[j] * timeStep;
            }
            //cout << cube[i].p[j] << endl;
        //gravityenergy
            gravityenergy += joint[i].m * gravity * joint[i].p[2];
            //kinetic
            double norm_v = sqrt(pow(joint[i].v[0], 2) + pow(joint[i].v[1], 2) + pow(joint[i].v[2], 2));
            kineticenergy += joint[i].m * pow(norm_v, 2) / 2;


            //            cout << i << joint[i].p[2] << endl;
            //            cout << cForces[1][2] << endl;
        }
        //cout << "hhhhh"<< endl;

        //cout << "wtfh" << endl;

        drawground();

        T = T + timeStep;

    }

};

vector<vector<GENE>> populationGene;
vector<Robot> robot;
vector<double> totald2;
vector<vector<GENE>> newgenepop;
vector<vector<GENE>> nextpop;
vector<GENE> bestGene;
int populationsize = size;



vector<double> bestdistance;

//double maxdistance;
vector<vector<GENE>> bestgene;


    vector<GENE> crossover(vector<GENE> p1, vector<GENE> p2) {
        if(p1.size()<p2.size()){
            int crossposition = 7 + rand() % static_cast<int>(p1.size());
                    
            vector<GENE> offspring = p1;
            //cout << "crossover" << endl;
            for (int i = 0; i < crossposition; i++) {
                offspring[i] = p2[i];
            }

            return offspring;
            }
        
        else{
            int crossposition = 7 + rand() % static_cast<int>(p2.size());
            vector<GENE> offspring = p2;
            //cout << "crossover" << endl;
            for (int i = 0; i < crossposition; i++) {
                offspring[i] = p1[i];
            }

            return offspring;
        }
    }

    vector<vector<GENE>> genereateGene(int num_foot) {
        vector<vector<GENE>> populationGene1;
        for (int i = 0; i < populationsize; i++) {
            vector<GENE> temp;
            GENE temp_2;
            double k_1 = 9000;
            double b_1 = 0;
            double c_1 = 0;
            int num_sp = 15+6*num_foot;
            for (int j = 0; j < num_sp; j++) {
                double k_1 = 100 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (1000 - 100)));
//                if (j == 18 or j == 37 or j == 17 or j == 27 or j == 29 or j == 39 or j == 46 or j == 48 or j == 54 or j == 56 or j == 62 or j == 64) {
                if (j < num_sp && j > 15) {
                    k_1 = 600 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (5000 - 600)));
                    b_1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 1);
                    c_1 = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));
                    }
                else{
                    k_1 = 9000;
                    b_1 = 0;
                    c_1 = 0;
                }
//

                temp_2 = { k_1, b_1, c_1 };
                temp.push_back(temp_2);
                //                tempVec.push_back(k1);
                //                tempVec.push_back(b1);
                //                tempVec.push_back(c1);

            }
            populationGene1.push_back(temp);
            
    }
        
        return populationGene1;
    }
    vector<vector<GENE>> selection(vector<vector<GENE>> populationGene,vector<double> totald) {
        vector<vector<GENE>> newgenepop;
        vector<int> index;
        populationGene.push_back()
        for (int i = 0; i < populationsize; i++) {
            double temp_min = totald[0];
            int pointer = 0;
            for (int j = i; j < populationsize; j++) {
                if (totald[j] > temp_min) {
                    temp_min = totald[j];
                    pointer = j;
                }

            }
            index.push_back(pointer);
        }

        newgenepop.clear();
        newgenepop.shrink_to_fit();
        for (int i = 0; i < index.size() / 2; i++) {
            newgenepop.push_back(populationGene[index[i]]);
        }
        bestgene.push_back(populationGene[index[0]]);
        bestdistance.push_back(totald[index[0]]);
        cout<<"best distance"<< totald[index[0]]<<endl;
        evaluation++;
//
//        outFile1<<evaluation<<" "<< totald[index[0]]<<endl;

        cout<<"evaluation: "<<evaluation<<endl;
//        cout<<"worst"<<totald[index[-1]]<<endl;
        //        for (int i=0;i<newgenepop.size(); i++){
        //            bestGene << newgenepop[i].k << " " <<newgenepop[i].b<< " " <<newgenepop[i].c << " ";
        //        }
        //        bestGene << "\n"
        return newgenepop;
    }
vector<GENE> mutation(vector<GENE> p1) {
    vector<int> list={};
    for(int i=15; i<p1.size();i++){
        list.push_back(i);
    }
    int i = 15 + rand() % static_cast<int>(p1.size()-15);
    int i1 = list[i];
    vector<GENE> offspring = p1;
    GENE temp;
    temp.k = 600 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (1000 - 100)));
    temp.b = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 1);
    temp.c = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));

    offspring[i1] = temp;

    return offspring;
}
    vector<vector<GENE>> nextgeneration(vector<vector<GENE>> newgenepop) {
        vector<vector<GENE>> nextpop = newgenepop;
        vector<GENE> p3;
        for (int a = 0; a<int(newgenepop.size()); a++) {
            int i = rand() % static_cast<int>(newgenepop.size());
            int j = rand() % static_cast<int>(newgenepop.size());

            vector<GENE> p1 = newgenepop[i];
            vector<GENE> p2 = newgenepop[j];
            if(p1.size() < p2.size()){
                p1 = crossover(p1, p2);
                nextpop.push_back(p1);
            }
            else{
                p2 = crossover(p1, p2);
                nextpop.push_back(p2);
            }
                
        }

        double r = static_cast <float> (rand()) / (static_cast <float> (RAND_MAX));
        if (r < 0.5) {
            vector<GENE> temp_off;
            int ra_ro = 1 + rand() % static_cast<int>(nextpop.size()-1);
            temp_off = mutation(nextpop[ra_ro]);
            nextpop[ra_ro] = temp_off;
            cout<<"mutationing"<<endl;


        }
//        for(int i =0; i<10; i++){
//            outFile2<<evaluation<<" "<<i <<" "<< nextpop[0][i].k<<" "<<nextpop[0][i].b<<" "<<nextpop[0][i].c<<endl;
//        }
        return nextpop;
    }

vector<Robot> getnewrobot(vector<Robot> robots, vector<vector<GENE>> nextpop) {
        
            //double x = -2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 4));
            double y = -1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2));
            //robots.push_back(ROBOT(0.0, 3.0*(i-populationSize/2), 0.0, populationGene[i]));
//            int num = 2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 8));
            int num_foot = 2+rand()%20;
            vector<GENE> temp;
                GENE temp_2;
                double k_1 = 9000;
                double b_1 = 0;
                double c_1 = 0;
                int num_sp = 15+6*num_foot;
                for (int j = 0; j < num_sp; j++) {
                    double k_1 = 100 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (1000 - 100)));
    //                if (j == 18 or j == 37 or j == 17 or j == 27 or j == 29 or j == 39 or j == 46 or j == 48 or j == 54 or j == 56 or j == 62 or j == 64) {
                    if (j < num_sp && j > 15) {
                        k_1 = 600 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (5000 - 600)));
                        b_1 = static_cast <float> (rand()) / static_cast <float> (RAND_MAX / 1);
                        c_1 = -2 * M_PI + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (4 * M_PI)));
                        }
                    else{
                        k_1 = 9000;
                        b_1 = 0;
                        c_1 = 0;
                    }
    //

                    temp_2 = { k_1, b_1, c_1 };
                    temp.push_back(temp_2);
            Robot robot1(0, y, 0.2, num_foot, temp);
            robots.push_back(robot1);

        
        return robots;
    }

    vector<Robot> getinitalrobot(vector<Robot> robots, vector<vector<GENE>> nextpop, int num_foot) {
        for (int i = 0; i < populationsize; i++) {

            //double x = -2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 4));
            double y = -1 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 2));
            //robots.push_back(ROBOT(0.0, 3.0*(i-populationSize/2), 0.0, populationGene[i]));
//            int num = 2 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / 8));
            
            Robot robot1(0, y, 0.2, num_foot, nextpop[i]);
            robots.push_back(robot1);

        }
        return robots;
    }
//    void move() {
//        for (int i = 0; i < populationsize; i++) {
//            robots[i].simulate();
//            //robots[i].drawplain();
//        }
//
//    }
//    void draw() {
//        for (int i = 0; i < populationsize; i++) {
//            robots[i].drawcube();
//
//        }
//
//    }
    double fittness(Robot r, int num_m) {
        double fit = 0;
        double startpoint;
        double sumdis = 0;
        for (int i = 0; i < num_m; i++) {
            fit += (r.joint[i].p[0] - r.initl[0]);
            sumdis += r.joint[i].p[0] - r.initl[0];
        }
        startpoint = sumdis/num_m;
        r.initl[0] = startpoint;
//        outFile3
//
//
//        << evaluation << " " << fit/13 << endl;
        //cout<<"startpoint"<<startpoint<<endl;
        return fit / num_m;
    }


    vector<double> totaldistance(vector<Robot> robots) {
        vector<double> totald;
        for (int i = 0; i < populationsize; i++) {
            double f = fittness(robots[i], robots[i].num_foot+7);
            totald.push_back(f);
            
        }

        return totald;
    }







void Print(const char* format, ...)
{
    char    buf[LEN];
    char* ch = buf;
    va_list args;
    //  Turn the parameters into a character string
    va_start(args, format);
    vsnprintf(buf, LEN, format, args);
    va_end(args);
    //  Display the characters one at a time at the current raster position
    while (*ch)
        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *ch++);
}


 *  OpenGL (GLUT) calls this routine to display the scene
 


void display()
{
    //drawground();
    double len = 0.2;  //  Length of axes
    //  Erase the window and the de5th buffer
    glClearColor(0.5372549, 0.6549019, 0.760784, 1.0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //  Enable Z-buffering in OpenGL
    glEnable(GL_DEPTH_TEST);
    //  Undo previous transformations
    glLoadIdentity();
    //  Eye position
    double Ex = -1 * dim * Sin(th) * Cos(ph);
    double Ey = +1 * dim * Sin(ph);
    double Ez = +1 * dim * Cos(th) * Cos(ph);
    gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, Cos(ph), 0);

    //Simulate();
    //loop
    //for (int i = 0; i < generation; i++) {
        if (T<6) {

            for (int j = 0; j<populationsize; j++){

            robot[j].simulate();
            
                robot[j].drawcube();}
            
        
         }
        else {
            totald2 = totaldistance(robot);
            newgenepop = selection(populationGene,totald2);
            nextpop = nextgeneration(newgenepop);
            populationGene = nextpop;
            //cout<<"in loop evolve"<<endl;
            T = 0;
        }


    //Draw axes
    //glColor3f(1, 0, 0);
//    if (axes)
//    {
//        glBegin(GL_LINES);
//        glLineWidth(2);
//        glVertex3d(0.0, 0.0, 0.0);
//        glVertex3d(len, 0.0, 0.0);
//        glVertex3d(0.0, 0.0, 0.0);
//        glVertex3d(0.0, len, 0.0);
//        glVertex3d(0.0, 0.0, 0.0);
//        glVertex3d(0.0, 0.0, len);
//        glEnd();
//        //  Label axes
//        glRasterPos3d(len, 0.0, 0.0);
//        Print("X");
//        glRasterPos3d(0.0, len, 0.0);
//        Print("Y");
//        glRasterPos3d(0.0, 0.0, len);
//        Print("Z");
//    }
    //  Render the scene
    glFlush();
    //  Make the rendered scene visible
    glutSwapBuffers();
}

/*
 *  GLUT calls this routine when an arrow key is pressed
 */
//void special(int key, int x, int y)
//{
//    //  Right arrow key - increase angle by 5 degrees
//    if (key == GLUT_KEY_RIGHT)
//        th += 5;
//    //  Left arrow key - decrease angle by 5 degrees
//    else if (key == GLUT_KEY_LEFT)
//        th -= 5;
//    //  Up arrow key - increase elevation by 5 degrees
//    else if (key == GLUT_KEY_UP)
//    {
//        if (ph + 5 < 90)
//        {
//            ph += 5;
//        }
//    }
//    //  Down arrow key - decrease elevation by 5 degrees
//    else if (key == GLUT_KEY_DOWN)
//    {
//        if (ph - 5 > 0)
//        {
//            ph -= 5;
//        }
//    }
//    //  Keep angles to +/-360 degrees
//    th %= 360;
//    ph %= 360;
//    //  Tell GLUT it is necessary to redisplay the scene
//    glutPostRedisplay();
//}

/*
 *  Set projection
 */
/*
void Project(double fov, double asp, double dim)
{
    //  Tell OpenGL we want to manipulate the projection matrix
    glMatrixMode(GL_PROJECTION);
    //  Undo previous transformations
    glLoadIdentity();
    //  Perspective transformation
    if (fov)
        gluPerspective(fov, asp, dim / 16, 16 * dim);
    //  Orthogonal transformation
    else
        glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);
    //  Switch to manipulating the model matrix
    glMatrixMode(GL_MODELVIEW);
    //  Undo previous transformations
    glLoadIdentity();
}

/*
 *  GLUT calls this routine when a key is pressed
 */
/*
void key(unsigned char ch, int x, int y)
{
    //  Exit on ESC
    if (ch == 27)
        exit(0);
    //  Reset view angle
    else if (ch == '0')
        th = ph = 0;
    //  Toggle axes
    else if (ch == 'a' || ch == 'A')
        axes = 1 - axes;
    //  Change field of view angle
    else if (ch == '-' && ch > 1)
        fov++;
    else if (ch == '=' && ch < 179)
        fov--;
    //  PageUp key - increase dim
    else if (ch == GLUT_KEY_PAGE_DOWN) {
        dim += 0.1;
    }
    //  PageDown key - decrease dim
    else if (ch == GLUT_KEY_PAGE_UP && dim > 1) {
        dim -= 0.1;
    }
    //  Keep angles to +/-360 degrees
    th %= 360;
    ph %= 360;
    //  Reproject
    Project(fov, asp, dim);
    //  Tell GLUT it is necessary to redisplay the scene
    glutPostRedisplay();
}

/*
 *  GLUT calls this routine when the window is resized
 */
/*
void reshape(int width, int height)
{
    //  Ratio of the width to the height of the window
    asp = (height > 0) ? (double)width / height : 1;
    //  Set the viewport to the entire window
    glViewport(0, 0, width, height);
    //  Set projection
    Project(fov, asp, dim);
}

/*
 *  GLUT calls this toutine when there is nothing else to do
 */
/*
void idle()
{
    glutPostRedisplay();
}

int main(int argc, char* argv[])
{
    
    
    int foot = 5;
    
    populationGene = genereateGene(foot);
    robot = getrobot(robot,populationGene);
    
    // Initialize GLUT and process user parameters
    glutInit(&argc, argv);
    // double buffered, true color 600*600
    glutInitWindowSize(1000, 800);
    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
    // create the window
    glutCreateWindow("evolve");
    //  Tell GLUT to call "idle" when there is nothing else to do
    glutIdleFunc(idle);
    //  Tell GLUT to call "display" when the scene should be drawn
    glutDisplayFunc(display);
    //  Tell GLUT to call "reshape" when the window is resized
    glutReshapeFunc(reshape);
    //  Tell GLUT to call "special" when an arrow key is pressed
    glutSpecialFunc(special);
    //  Tell GLUT to call "key" when a key is pressed
    glutKeyboardFunc(key);
    init();

    //  Pass control to GLUT so it can interact with the user
    glutMainLoop();
        
    return 0;
    
    

};
*/

//int ppp = 0;
//GLfloat worldRotation[16] = { 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1 };
///*std::clock_t start = std::clock();
//double duration;
//
//std::ofstream myFile("breathing.txt");*/
//ofstream outfile4("efficienty.txt");
//struct MASS
//{
//    double m;       // mass
//    double p[3];    // 3D position
//    double v[3];    // 3D velocity
//    double a[3];    // 3D acceleration
//};
//
//struct SPRING
//{
//    double k;       // spring constant
//    double L_0;     // rest length
//    int m1;         // first mass connected
//    int m2;         // second mass connected
//};
//
//vector<MASS> cubemass(double mass, double length, double x, double y, double z)
//{
//    vector<MASS> cube(8);
//    cube[0] = { mass, {x,y,z}, {0,0,0}, {0,0,0} };
//    cube[1] = { mass, {x ,y+length,z}, {0,0,0}, {0,0,0} };
//    cube[2] = { mass, {x ,y,z+length}, {0,0,0}, {0,0,0} };
//    cube[3] = { mass, {x + length,y + length,z}, {0,0,0}, {0,0,0} };
//    cube[4] = { mass, {x + length,y,z + length}, {0,0,0}, {0,0,0} };
//    cube[5] = { mass, {x,y + length,z + length}, {0,0,0}, {0,0,0} };
//    cube[6] = { mass, {x + length,y,z}, {0,0,0}, {0,0,0} };
//    cube[7] = { mass, {x + length,y + length,z + length}, {0,0,0}, {0,0,0} };
//    return cube;
//}
//vector<SPRING> cubespring(double length, double k)
//{
//    double short_diagonals = sqrt(2) * length;
//    double long_diagonals = sqrt(3) * length;
//    vector<SPRING> spring(28);
//    spring[0] = { k,length,0,1 };
//    spring[1] = { k,length,0,6 };
//    spring[2] = { k,length,0,2 };
//    spring[3] = { k,length,1,3 };
//    spring[4] = { k,length,3,6 };
//    spring[5] = { k,length,2,4 };
//    spring[6] = { k,length,4,6 };
//    spring[7] = { k,length,2,5};
//    spring[8] = { k,length,1,5 };
//    spring[9] = { k,length,5,7 };
//    spring[10] = { k,length,4,7};
//    spring[11] = { k,length,3,7 };
//
//    spring[12] = { k,short_diagonals,3,4 };
//    spring[13] = { k,short_diagonals,6,7 };
//    spring[14] = { k,short_diagonals,1,7 };
//    spring[15] = { k,short_diagonals,3,5 };
//    spring[16] = { k,short_diagonals,1,2 };
//    spring[17] = { k,short_diagonals,0,5 };
//    spring[18] = { k,short_diagonals,2,6 };
//    spring[19] = { k,short_diagonals,0,4 };
//    spring[20] = { k,short_diagonals,2,7 };
//    spring[21] = { k,short_diagonals,4,5 };
//    spring[22] = { k,short_diagonals,0,3 };
//    spring[23] = { k,short_diagonals,1,6 };
//
//    spring[24] = { k,long_diagonals,0,7 };
//    spring[25] = { k,long_diagonals,2,3 };
//    spring[26] = { k,long_diagonals,1,4 };
//    spring[27] = { k,long_diagonals,5,6 };
//    return spring;
//}
//vector<MASS> cube = cubemass(mass, length, 0, 0, 0);
//vector<vector<double>> cForces(8, vector<double>(3));
//vector<SPRING> spring = cubespring(length, 1000);
//GLuint tex;
//GLUquadric* sphere;
//
//void make_tex(void)
//{
//    unsigned char data[256][256][3];
//    for (int y = 0; y < 255; y++) {
//        for (int x = 0; x < 255; x++) {
//            unsigned char* p = data[y][x];
//            p[0] = p[1] = p[2] = (x ^ y) & 8 ? 255 : 0;
//        }
//    }
//    glGenTextures(1, &tex);
//    glBindTexture(GL_TEXTURE_2D, tex);
//    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, (const GLvoid*)data);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//}
//
//void init(void)
//{
//    glEnable(GL_DEPTH_TEST);
//    make_tex();
//    sphere = gluNewQuadric();
//    glEnable(GL_TEXTURE_2D);
//}
//
//void drawground()
//{
//    glPushMatrix();
//    glColor3f(0.96078,0.96078,0.86274);
//    glBindTexture(GL_TEXTURE_2D,grassTexture);
//    glBegin(GL_QUADS);
//    glNormal3f( 0, 1, 0);
//    glTexCoord2f(0.0,0.0);  glVertex3f(-0.5,+0.0,-0.5);
//    glTexCoord2f(0.0,1.0);  glVertex3f(+0.5,+0.0,-0.5);
//    glTexCoord2f(1.0,1.0);  glVertex3f(+0.5,+0.0,+0.5);
//    glTexCoord2f(1.0,0.0);  glVertex3f(-0.5,+0.0,+0.5);
//    glEnd();
//    glPopMatrix();
//    glDisable(GL_TEXTURE_2D);
//    for (int i=0; i < 10; i++) {
//   for (int j=-9;  j<10; j++) {
//        glColor3f(0, 1, 0);
//        glPushMatrix();
//        glMultMatrixf(worldRotation);
//        glBegin(GL_LINES);
//        glLineWidth(10);
//        glVertex3f(-0.5*i/10,0.02 , 0.5*j/10);
////        glV/Users/chiqu/class2019fall/EA/assignment3/HW3/HW3.cppertex3f(0.5*i/10, 0.5*j/10, 0.01);
//        glEnd();
//        glPopMatrix();
//   }}
//}
//void drawcube()
//{
//
//
////    GLUquadric* quad[8];
////
////    for (int i=0; i < 8; i++) {
////        glColor3f(0.7, 0.3, 0.3);
////        quad[i] = gluNewQuadric();
////        glPushMatrix();
////        glMultMatrixf(worldRotation);
////        glTranslated(cube[i].p[0], cube[i].p[1], cube[i].p[2]);
////
////        gluQuadricDrawStyle(quad[i], GLU_FILL);
////        glBindTexture(GL_TEXTURE_2D, tex);
////        gluQuadricTexture(quad[i], GL_TRUE);
////        gluQuadricNormals(quad[i], GLU_SMOOTH);
////        gluSphere(quad[i], 1.0 / 100, 10, 10);
////        glPopMatrix();
////    }
//
//    for (int i=0; i < 8; i++) {
//        for (int j=i+1;  j<8; j++) {
//            glColor3f(0.2*i, 0.3*i, 0.5*i);
//            glPushMatrix();
//            glMultMatrixf(worldRotation);
//            glBegin(GL_LINES);
//            glLineWidth(10);
//            glVertex3f(cube[i].p[0], cube[i].p[1], cube[i].p[2]);
//            glVertex3f(cube[j].p[0], cube[j].p[1], cube[j].p[2]);
//            glEnd();
//            glPopMatrix();
//
//        }
//    }
//
//       glPushMatrix();
//       glMultMatrixf(worldRotation);
//
//
//
//       //  Front
//       glColor3f(242.0/255, 138.0/255, 58.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[4].p[0],cube[4].p[1],cube[4].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[7].p[0],cube[7].p[1],cube[7].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[3].p[0],cube[3].p[1],cube[3].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[6].p[0],cube[6].p[1],cube[6].p[2]);
//       glEnd();
//
//
//    //  Back
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[2].p[0],cube[2].p[1],cube[2].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[5].p[0],cube[5].p[1],cube[5].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[1].p[0],cube[1].p[1],cube[1].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[0].p[0],cube[0].p[1],cube[0].p[2]);
//       glEnd();
//
//    //  Right
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[5].p[0],cube[5].p[1],cube[5].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[7].p[0],cube[7].p[1],cube[7].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[3].p[0],cube[3].p[1],cube[3].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[1].p[0],cube[1].p[1],cube[1].p[2]);
//       glEnd();
//
//    //  Left
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[2].p[0],cube[2].p[1],cube[2].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[4].p[0],cube[4].p[1],cube[4].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[6].p[0],cube[6].p[1],cube[6].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[0].p[0],cube[0].p[1],cube[0].p[2]);
//       glEnd();
//
//    //  Top
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[5].p[0],cube[5].p[1],cube[5].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[7].p[0],cube[7].p[1],cube[7].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[4].p[0],cube[4].p[1],cube[4].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[2].p[0],cube[2].p[1],cube[2].p[2]);
//       glEnd();
//
//    //  Bottom
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[0].p[0],cube[0].p[1],cube[0].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[1].p[0],cube[1].p[1],cube[1].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[3].p[0],cube[3].p[1],cube[3].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[6].p[0],cube[6].p[1],cube[6].p[2]);
//       glEnd();
//        glPopMatrix();
//     glDisable(GL_TEXTURE_2D);
//
//
//
//
//}
//
//float L(MASS mass1, MASS mass2) {
//    double length = sqrt(pow((mass1.p[0] - mass2.p[0]), 2) + pow((mass1.p[1] - mass2.p[1]), 2) + pow((mass1.p[2] - mass2.p[2]), 2));
//
//    return length;
//}
//void simulate() {
//
//    for (int i = 0; i < 8; i++) {
//    cForces[i][0] = 0.0;
//    cForces[i][1] = 0.0;
//    cForces[i][2] = -cube[i].m * gravity;}
////    if (oneforce==0){
////        cForces[1][0] = 0.1;
////        oneforce = 1;
////    }
//
//    if (T == 0.005) {
//     cForces[5][1] = 2.0;
//    }
//
//    //cout << "hhhhh" <<cForces[1][2]<< endl;
//    //cout << spring[2].L_0 << endl;
//    for (int i = 0; i < 28; i++) {
//        cout<<ppp<<endl;
//
//        ppp++;
//
//        if (T==0.001 or T == 1.0 or T == 2.0 or T == 3.0) {
//         outfile4 << " " << T << " " << ppp << endl;
//         ppp = 0;
//        }
//        cout<<T<<endl;
//        if (T > 0) {
//                spring[24].L_0 = 1.0 * length + 0.05 * length * sin(20 * T);
//                    spring[25].L_0 = 1.0 * length + 0.05 * length * sin(20 * T);
//                   spring[26].L_0 = 1.0 * length + 0.08 * length * sin(20 * T);
//                   spring[27].L_0 = 1.0 * length + 0.08 * length * sin(20 * T);
//        }
//        MASS mass1 = cube[spring[i].m1];
//        MASS mass2 = cube[spring[i].m2];
//        double pd[3] = { mass2.p[0] - mass1.p[0],mass2.p[1] - mass1.p[1],mass2.p[2] - mass1.p[2] };
//        double new_L = L( mass1, mass2);
//        double L_0 = spring[i].L_0;
//        double force = k * fabs(new_L - L_0);
//        //cout <<i<<"---new_L---" <<new_L << endl;
//        double norm_pd[3] = { pd[0] / new_L, pd[1] / new_L, pd[2] / new_L };
//        //cout << i << endl;
//        //
//
//        //cout << i<<"---force---" <<force << endl;
//        //compression
//        if (new_L < spring[i].L_0) {
//            cForces[spring[i].m1][0] -= norm_pd[0] * force;
//            cForces[spring[i].m1][1] -= norm_pd[1] * force;
//            cForces[spring[i].m1][2] -= norm_pd[2] * force;
//            cForces[spring[i].m2][0] += norm_pd[0] * force;
//            cForces[spring[i].m2][1] += norm_pd[1] * force;
//            cForces[spring[i].m2][2] += norm_pd[2] * force;
//        }
//
//        //tension
//        else{
//            cForces[spring[i].m1][0] += norm_pd[0] * force;
//            cForces[spring[i].m1][1] += norm_pd[1] * force;
//            cForces[spring[i].m1][2] += norm_pd[2] * force;
//            cForces[spring[i].m2][0] -= norm_pd[0] * force;
//            cForces[spring[i].m2][1] -= norm_pd[1] * force;
//            cForces[spring[i].m2][2] -= norm_pd[2] * force;
//        }
//        //cout <<i<<"---52---" <<cForces[5][2] << endl;
//    }
//    //cout << "force" << endl;
//    //cout << cForces[2][2] << endl;
//    for (int i = 0; i < 8; i++) {
//        //cout << cForces[i][2] << endl;
//        //cout << 's' << endl;
//        if (cube[i].p[2] < 0) {
//         cForces[i][2] -= Nground * cube[i].p[2];
//         double Fh = sqrt(pow(cForces[i][0], 2) + pow(cForces[i][1], 2));
//         double Fv = cForces[i][2];
//         if (Fh < Fv * frictionCoefficient) {
//          cForces[i][0] = 0;
//          cForces[i][1] = 0;
//          cube[i].v[0] = 0;
//          cube[i].v[1] = 0;
//         }
//         else {
//          double Fh_new = Fh -  Fv * frictionCoefficient;
//          cForces[i][0] = cForces[i][0] * Fh_new / Fh;
//          cForces[i][1] = cForces[i][1] * Fh_new / Fh;
//         }
//        }
////        if (cube[i].p[2] < 0) {
////            cForces[i][2] -= Nground * cube[i].p[2];
////        }
//        for (int j = 0; j < 3; j++) {
//            cube[i].a[j] = cForces[i][j] / cube[i].m;
//            cube[i].v[j] += cube[i].a[j] * timeStep;
//            cube[i].p[j] += cube[i].v[j] * timeStep;
//            //cout << cube[i].p[j] << endl;
//
//        }
//    //cout <<i <<cube[i].p[2] << endl;
//    //cout << cForces[1][2] << endl;
//    }
//    drawcube();
//    drawground();
//    T = T + timeStep;
//}
//
//
//
//void Print(const char* format, ...)
//{
//    char    buf[LEN];
//    char* ch = buf;
//    va_list args;
//    //  Turn the parameters into a character string
//    va_start(args, format);
//    vsnprintf(buf, LEN, format, args);
//    va_end(args);
//    //  Display the characters one at a time at the current raster position
//    while (*ch)
//        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *ch++);
//}
//
///*
// *  OpenGL (GLUT) calls this routine to display the scene
// */
//void display()
//{
//
//    const double len = 0.2;  //  Length of axes
//    //  Erase the window and the depth buffer
//    glClearColor(0.5372549, 0.6549019, 0.760784, 1.0);
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    //  Enable Z-buffering in OpenGL
//    glEnable(GL_DEPTH_TEST);
//    //  Undo previous transformations
//    glLoadIdentity();
//    //  Eye position
//    double Ex = -1 * dim * Sin(th) * Cos(ph);
//    double Ey = +1 * dim * Sin(ph);
//    double Ez = +1 * dim * Cos(th) * Cos(ph);
//    gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, Cos(ph), 0);
//
//    simulate();
//    //drawcube();
//
//
//    //  Draw axes
//    glColor3f(1, 0, 0);
////    if (axes)
////    {
////        glBegin(GL_LINES);
////        glLineWidth(2);
////        glVertex3d(0.0, 0.0, 0.0);
////        glVertex3d(len, 0.0, 0.0);
////        glVertex3d(0.0, 0.0, 0.0);
////        glVertex3d(0.0, len, 0.0);
////        glVertex3d(0.0, 0.0, 0.0);
////        glVertex3d(0.0, 0.0, len);
////        glEnd();
////        //  Label axes
////        glRasterPos3d(len, 0.0, 0.0);
////        Print("X");
////        glRasterPos3d(0.0, len, 0.0);
////        Print("Y");
////        glRasterPos3d(0.0, 0.0, len);
////        Print("Z");
////    }
//    //  Render the scene
//    glFlush();
//    //  Make the rendered scene visible
//    glutSwapBuffers();
//}
//
///*
// *  GLUT calls this routine when an arrow key is pressed
// */
//void special(int key, int x, int y)
//{
//    //  Right arrow key - increase angle by 5 degrees
//    if (key == GLUT_KEY_RIGHT)
//        th += 5;
//    //  Left arrow key - decrease angle by 5 degrees
//    else if (key == GLUT_KEY_LEFT)
//        th -= 5;
//    //  Up arrow key - increase elevation by 5 degrees
//    else if (key == GLUT_KEY_UP)
//    {
//        if (ph + 5 < 90)
//        {
//            ph += 5;
//        }
//    }
//    //  Down arrow key - decrease elevation by 5 degrees
//    else if (key == GLUT_KEY_DOWN)
//    {
//        if (ph - 5 > 0)
//        {
//            ph -= 5;
//        }
//    }
//    //  Keep angles to +/-360 degrees
//    th %= 360;
//    ph %= 360;
//    //  Tell GLUT it is necessary to redisplay the scene
//    glutPostRedisplay();
//}
//
///*
// *  Set projection
// */
//void Project(double fov, double asp, double dim)
//{
//    //  Tell OpenGL we want to manipulate the projection matrix
//    glMatrixMode(GL_PROJECTION);
//    //  Undo previous transformations
//    glLoadIdentity();
//    //  Perspective transformation
//    if (fov)
//        gluPerspective(fov, asp, dim / 16, 16 * dim);
//    //  Orthogonal transformation
//    else
//        glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);
//    //  Switch to manipulating the model matrix
//    glMatrixMode(GL_MODELVIEW);
//    //  Undo previous transformations
//    glLoadIdentity();
//}
//
///*
// *  GLUT calls this routine when a key is pressed
// */
//void key(unsigned char ch, int x, int y)
//{
//    //  Exit on ESC
//    if (ch == 27)
//        exit(0);
//    //  Reset view angle
//    else if (ch == '0')
//        th = ph = 0;
//    //  Toggle axes
//    else if (ch == 'a' || ch == 'A')
//        axes = 1 - axes;
//    //  Change field of view angle
//    else if (ch == '-' && ch > 1)
//        fov++;
//    else if (ch == '=' && ch < 179)
//        fov--;
//    //  PageUp key - increase dim
//    else if (ch == GLUT_KEY_PAGE_DOWN) {
//        dim += 0.1;
//    }
//    //  PageDown key - decrease dim
//    else if (ch == GLUT_KEY_PAGE_UP && dim > 1) {
//        dim -= 0.1;
//    }
//    //  Keep angles to +/-360 degrees
//    th %= 360;
//    ph %= 360;
//    //  Reproject
//    Project(fov, asp, dim);
//    //  Tell GLUT it is necessary to redisplay the scene
//    glutPostRedisplay();
//}
//
///*
// *  GLUT calls this routine when the window is resized
// */
//void reshape(int width, int height)
//{
//    //  Ratio of the width to the height of the window
//    asp = (height > 0) ? (double)width / height : 1;
//    //  Set the viewport to the entire window
//    glViewport(0, 0, width, height);
//    //  Set projection
//    Project(fov, asp, dim);
//}
//
///*
// *  GLUT calls this toutine when there is nothing else to do
// */
//void idle()
//{
//    glutPostRedisplay();
//}
//
//int main(int argc, char* argv[])
//{
//    // Initialize GLUT and process user parameters
//    glutInit(&argc, argv);
//    // double buffered, true color 600*600
//    glutInitWindowSize(1000, 800);
//    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
//    // create the window
//    glutCreateWindow("Slight Spin_yj2563_cl3895");
//    //  Tell GLUT to call "idle" when there is nothing else to do
//    glutIdleFunc(idle);
//    //  Tell GLUT to call "display" when the scene should be drawn
//    glutDisplayFunc(display);
//    //  Tell GLUT to call "reshape" when the window is resized
//    glutReshapeFunc(reshape);
//    //  Tell GLUT to call "special" when an arrow key is pressed
//    glutSpecialFunc(special);
//    //  Tell GLUT to call "key" when a key is pressed
//    glutKeyboardFunc(key);
//    init();
//    //  Pass control to GLUT so it can interact with the user
//    glutMainLoop();
//    return 0;
//};
//

//int ppp = 0;
//GLfloat worldRotation[16] = { 1,0,0,0,0,0,1,0,0,1,0,0,0,0,0,1 };
///*std::clock_t start = std::clock();
//double duration;
//
//std::ofstream myFile("breathing.txt");*/
//ofstream outfile4("efficienty.txt");
//struct MASS
//{
//    double m;       // mass
//    double p[3];    // 3D position
//    double v[3];    // 3D velocity
//    double a[3];    // 3D acceleration
//};
//
//struct SPRING
//{
//    double k;       // spring constant
//    double L_0;     // rest length
//    int m1;         // first mass connected
//    int m2;         // second mass connected
//};
//
//vector<MASS> cubemass(double mass, double length, double x, double y, double z)
//{
//    vector<MASS> cube(8);
//    cube[0] = { mass, {x,y,z}, {0,0,0}, {0,0,0} };
//    cube[1] = { mass, {x ,y+length,z}, {0,0,0}, {0,0,0} };
//    cube[2] = { mass, {x ,y,z+length}, {0,0,0}, {0,0,0} };
//    cube[3] = { mass, {x + length,y + length,z}, {0,0,0}, {0,0,0} };
//    cube[4] = { mass, {x + length,y,z + length}, {0,0,0}, {0,0,0} };
//    cube[5] = { mass, {x,y + length,z + length}, {0,0,0}, {0,0,0} };
//    cube[6] = { mass, {x + length,y,z}, {0,0,0}, {0,0,0} };
//    cube[7] = { mass, {x + length,y + length,z + length}, {0,0,0}, {0,0,0} };
//    return cube;
//}
//vector<SPRING> cubespring(double length, double k)
//{
//    double short_diagonals = sqrt(2) * length;
//    double long_diagonals = sqrt(3) * length;
//    vector<SPRING> spring(28);
//    spring[0] = { k,length,0,1 };
//    spring[1] = { k,length,0,6 };
//    spring[2] = { k,length,0,2 };
//    spring[3] = { k,length,1,3 };
//    spring[4] = { k,length,3,6 };
//    spring[5] = { k,length,2,4 };
//    spring[6] = { k,length,4,6 };
//    spring[7] = { k,length,2,5};
//    spring[8] = { k,length,1,5 };
//    spring[9] = { k,length,5,7 };
//    spring[10] = { k,length,4,7};
//    spring[11] = { k,length,3,7 };
//
//    spring[12] = { k,short_diagonals,3,4 };
//    spring[13] = { k,short_diagonals,6,7 };
//    spring[14] = { k,short_diagonals,1,7 };
//    spring[15] = { k,short_diagonals,3,5 };
//    spring[16] = { k,short_diagonals,1,2 };
//    spring[17] = { k,short_diagonals,0,5 };
//    spring[18] = { k,short_diagonals,2,6 };
//    spring[19] = { k,short_diagonals,0,4 };
//    spring[20] = { k,short_diagonals,2,7 };
//    spring[21] = { k,short_diagonals,4,5 };
//    spring[22] = { k,short_diagonals,0,3 };
//    spring[23] = { k,short_diagonals,1,6 };
//
//    spring[24] = { k,long_diagonals,0,7 };
//    spring[25] = { k,long_diagonals,2,3 };
//    spring[26] = { k,long_diagonals,1,4 };
//    spring[27] = { k,long_diagonals,5,6 };
//    return spring;
//}
//vector<MASS> cube = cubemass(mass, length, 0, 0, 0);
//vector<vector<double>> cForces(8, vector<double>(3));
//vector<SPRING> spring = cubespring(length, 1000);
//GLuint tex;
//GLUquadric* sphere;
//
//void make_tex(void)
//{
//    unsigned char data[256][256][3];
//    for (int y = 0; y < 255; y++) {
//        for (int x = 0; x < 255; x++) {
//            unsigned char* p = data[y][x];
//            p[0] = p[1] = p[2] = (x ^ y) & 8 ? 255 : 0;
//        }
//    }
//    glGenTextures(1, &tex);
//    glBindTexture(GL_TEXTURE_2D, tex);
//    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 256, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, (const GLvoid*)data);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
//    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
//}
//
//void init(void)
//{
//    glEnable(GL_DEPTH_TEST);
//    make_tex();
//    sphere = gluNewQuadric();
//    glEnable(GL_TEXTURE_2D);
//}
//
//void drawground()
//{
//    glPushMatrix();
//    glColor3f(0.96078,0.96078,0.86274);
//    glBindTexture(GL_TEXTURE_2D,grassTexture);
//    glBegin(GL_QUADS);
//    glNormal3f( 0, 1, 0);
//    glTexCoord2f(0.0,0.0);  glVertex3f(-0.5,+0.0,-0.5);
//    glTexCoord2f(0.0,1.0);  glVertex3f(+0.5,+0.0,-0.5);
//    glTexCoord2f(1.0,1.0);  glVertex3f(+0.5,+0.0,+0.5);
//    glTexCoord2f(1.0,0.0);  glVertex3f(-0.5,+0.0,+0.5);
//    glEnd();
//    glPopMatrix();
//    glDisable(GL_TEXTURE_2D);
//    for (int i=0; i < 10; i++) {
//   for (int j=-9;  j<10; j++) {
//        glColor3f(0, 1, 0);
//        glPushMatrix();
//        glMultMatrixf(worldRotation);
//        glBegin(GL_LINES);
//        glLineWidth(10);
//        glVertex3f(-0.5*i/10,0.02 , 0.5*j/10);
////        glV/Users/chiqu/class2019fall/EA/assignment3/HW3/HW3.cppertex3f(0.5*i/10, 0.5*j/10, 0.01);
//        glEnd();
//        glPopMatrix();
//   }}
//}
//void drawcube()
//{
//
//
////    GLUquadric* quad[8];
////
////    for (int i=0; i < 8; i++) {
////        glColor3f(0.7, 0.3, 0.3);
////        quad[i] = gluNewQuadric();
////        glPushMatrix();
////        glMultMatrixf(worldRotation);
////        glTranslated(cube[i].p[0], cube[i].p[1], cube[i].p[2]);
////
////        gluQuadricDrawStyle(quad[i], GLU_FILL);
////        glBindTexture(GL_TEXTURE_2D, tex);
////        gluQuadricTexture(quad[i], GL_TRUE);
////        gluQuadricNormals(quad[i], GLU_SMOOTH);
////        gluSphere(quad[i], 1.0 / 100, 10, 10);
////        glPopMatrix();
////    }
//
//    for (int i=0; i < 8; i++) {
//        for (int j=i+1;  j<8; j++) {
//            glColor3f(0.2*i, 0.3*i, 0.5*i);
//            glPushMatrix();
//            glMultMatrixf(worldRotation);
//            glBegin(GL_LINES);
//            glLineWidth(10);
//            glVertex3f(cube[i].p[0], cube[i].p[1], cube[i].p[2]);
//            glVertex3f(cube[j].p[0], cube[j].p[1], cube[j].p[2]);
//            glEnd();
//            glPopMatrix();
//
//        }
//    }
//
//       glPushMatrix();
//       glMultMatrixf(worldRotation);
//
//
//
//       //  Front
//       glColor3f(242.0/255, 138.0/255, 58.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[4].p[0],cube[4].p[1],cube[4].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[7].p[0],cube[7].p[1],cube[7].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[3].p[0],cube[3].p[1],cube[3].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[6].p[0],cube[6].p[1],cube[6].p[2]);
//       glEnd();
//
//
//    //  Back
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[2].p[0],cube[2].p[1],cube[2].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[5].p[0],cube[5].p[1],cube[5].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[1].p[0],cube[1].p[1],cube[1].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[0].p[0],cube[0].p[1],cube[0].p[2]);
//       glEnd();
//
//    //  Right
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[5].p[0],cube[5].p[1],cube[5].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[7].p[0],cube[7].p[1],cube[7].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[3].p[0],cube[3].p[1],cube[3].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[1].p[0],cube[1].p[1],cube[1].p[2]);
//       glEnd();
//
//    //  Left
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[2].p[0],cube[2].p[1],cube[2].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[4].p[0],cube[4].p[1],cube[4].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[6].p[0],cube[6].p[1],cube[6].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[0].p[0],cube[0].p[1],cube[0].p[2]);
//       glEnd();
//
//    //  Top
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[5].p[0],cube[5].p[1],cube[5].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[7].p[0],cube[7].p[1],cube[7].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[4].p[0],cube[4].p[1],cube[4].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[2].p[0],cube[2].p[1],cube[2].p[2]);
//       glEnd();
//
//    //  Bottom
//       glColor3f(240.0/255, 136.0/255, 56.0/255);
//       glBindTexture(GL_TEXTURE_2D,slimeTexture);
//       glBegin(GL_QUADS);
//       glNormal3f( 0, 0, 1);
//       glTexCoord2f(0.0f,0.0f);    glVertex3f(cube[0].p[0],cube[0].p[1],cube[0].p[2]);
//       glTexCoord2f(1.0f,0.0f);    glVertex3f(cube[1].p[0],cube[1].p[1],cube[1].p[2]);
//       glTexCoord2f(1.0f,1.0f);    glVertex3f(cube[3].p[0],cube[3].p[1],cube[3].p[2]);
//       glTexCoord2f(0.0f,1.0f);    glVertex3f(cube[6].p[0],cube[6].p[1],cube[6].p[2]);
//       glEnd();
//        glPopMatrix();
//     glDisable(GL_TEXTURE_2D);
//
//
//
//
//}
//
//float L(MASS mass1, MASS mass2) {
//    double length = sqrt(pow((mass1.p[0] - mass2.p[0]), 2) + pow((mass1.p[1] - mass2.p[1]), 2) + pow((mass1.p[2] - mass2.p[2]), 2));
//
//    return length;
//}
//void simulate() {
//
//    for (int i = 0; i < 8; i++) {
//    cForces[i][0] = 0.0;
//    cForces[i][1] = 0.0;
//    cForces[i][2] = -cube[i].m * gravity;}
////    if (oneforce==0){
////        cForces[1][0] = 0.1;
////        oneforce = 1;
////    }
//
//    if (T == 0.005) {
//     cForces[5][1] = 2.0;
//    }
//
//    //cout << "hhhhh" <<cForces[1][2]<< endl;
//    //cout << spring[2].L_0 << endl;
//    for (int i = 0; i < 28; i++) {
//        cout<<ppp<<endl;
//
//        ppp++;
//
//        if (T==0.001 or T == 1.0 or T == 2.0 or T == 3.0) {
//         outfile4 << " " << T << " " << ppp << endl;
//         ppp = 0;
//        }
//        cout<<T<<endl;
//        if (T > 0) {
//                spring[24].L_0 = 1.0 * length + 0.05 * length * sin(20 * T);
//                    spring[25].L_0 = 1.0 * length + 0.05 * length * sin(20 * T);
//                   spring[26].L_0 = 1.0 * length + 0.08 * length * sin(20 * T);
//                   spring[27].L_0 = 1.0 * length + 0.08 * length * sin(20 * T);
//        }
//        MASS mass1 = cube[spring[i].m1];
//        MASS mass2 = cube[spring[i].m2];
//        double pd[3] = { mass2.p[0] - mass1.p[0],mass2.p[1] - mass1.p[1],mass2.p[2] - mass1.p[2] };
//        double new_L = L( mass1, mass2);
//        double L_0 = spring[i].L_0;
//        double force = k * fabs(new_L - L_0);
//        //cout <<i<<"---new_L---" <<new_L << endl;
//        double norm_pd[3] = { pd[0] / new_L, pd[1] / new_L, pd[2] / new_L };
//        //cout << i << endl;
//        //
//
//        //cout << i<<"---force---" <<force << endl;
//        //compression
//        if (new_L < spring[i].L_0) {
//            cForces[spring[i].m1][0] -= norm_pd[0] * force;
//            cForces[spring[i].m1][1] -= norm_pd[1] * force;
//            cForces[spring[i].m1][2] -= norm_pd[2] * force;
//            cForces[spring[i].m2][0] += norm_pd[0] * force;
//            cForces[spring[i].m2][1] += norm_pd[1] * force;
//            cForces[spring[i].m2][2] += norm_pd[2] * force;
//        }
//
//        //tension
//        else{
//            cForces[spring[i].m1][0] += norm_pd[0] * force;
//            cForces[spring[i].m1][1] += norm_pd[1] * force;
//            cForces[spring[i].m1][2] += norm_pd[2] * force;
//            cForces[spring[i].m2][0] -= norm_pd[0] * force;
//            cForces[spring[i].m2][1] -= norm_pd[1] * force;
//            cForces[spring[i].m2][2] -= norm_pd[2] * force;
//        }
//        //cout <<i<<"---52---" <<cForces[5][2] << endl;
//    }
//    //cout << "force" << endl;
//    //cout << cForces[2][2] << endl;
//    for (int i = 0; i < 8; i++) {
//        //cout << cForces[i][2] << endl;
//        //cout << 's' << endl;
//        if (cube[i].p[2] < 0) {
//         cForces[i][2] -= Nground * cube[i].p[2];
//         double Fh = sqrt(pow(cForces[i][0], 2) + pow(cForces[i][1], 2));
//         double Fv = cForces[i][2];
//         if (Fh < Fv * frictionCoefficient) {
//          cForces[i][0] = 0;
//          cForces[i][1] = 0;
//          cube[i].v[0] = 0;
//          cube[i].v[1] = 0;
//         }
//         else {
//          double Fh_new = Fh -  Fv * frictionCoefficient;
//          cForces[i][0] = cForces[i][0] * Fh_new / Fh;
//          cForces[i][1] = cForces[i][1] * Fh_new / Fh;
//         }
//        }
////        if (cube[i].p[2] < 0) {
////            cForces[i][2] -= Nground * cube[i].p[2];
////        }
//        for (int j = 0; j < 3; j++) {
//            cube[i].a[j] = cForces[i][j] / cube[i].m;
//            cube[i].v[j] += cube[i].a[j] * timeStep;
//            cube[i].p[j] += cube[i].v[j] * timeStep;
//            //cout << cube[i].p[j] << endl;
//
//        }
//    //cout <<i <<cube[i].p[2] << endl;
//    //cout << cForces[1][2] << endl;
//    }
//    drawcube();
//    drawground();
//    T = T + timeStep;
//}
//
//
//
//void Print(const char* format, ...)
//{
//    char    buf[LEN];
//    char* ch = buf;
//    va_list args;
//    //  Turn the parameters into a character string
//    va_start(args, format);
//    vsnprintf(buf, LEN, format, args);
//    va_end(args);
//    //  Display the characters one at a time at the current raster position
//    while (*ch)
//        glutBitmapCharacter(GLUT_BITMAP_HELVETICA_18, *ch++);
//}
//
///*
// *  OpenGL (GLUT) calls this routine to display the scene
// */
//void display()
//{
//
//    const double len = 0.2;  //  Length of axes
//    //  Erase the window and the depth buffer
//    glClearColor(0.5372549, 0.6549019, 0.760784, 1.0);
//    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
//    //  Enable Z-buffering in OpenGL
//    glEnable(GL_DEPTH_TEST);
//    //  Undo previous transformations
//    glLoadIdentity();
//    //  Eye position
//    double Ex = -1 * dim * Sin(th) * Cos(ph);
//    double Ey = +1 * dim * Sin(ph);
//    double Ez = +1 * dim * Cos(th) * Cos(ph);
//    gluLookAt(Ex, Ey, Ez, 0, 0, 0, 0, Cos(ph), 0);
//
//    simulate();
//    //drawcube();
//
//
//    //  Draw axes
//    glColor3f(1, 0, 0);
////    if (axes)
////    {
////        glBegin(GL_LINES);
////        glLineWidth(2);
////        glVertex3d(0.0, 0.0, 0.0);
////        glVertex3d(len, 0.0, 0.0);
////        glVertex3d(0.0, 0.0, 0.0);
////        glVertex3d(0.0, len, 0.0);
////        glVertex3d(0.0, 0.0, 0.0);
////        glVertex3d(0.0, 0.0, len);
////        glEnd();
////        //  Label axes
////        glRasterPos3d(len, 0.0, 0.0);
////        Print("X");
////        glRasterPos3d(0.0, len, 0.0);
////        Print("Y");
////        glRasterPos3d(0.0, 0.0, len);
////        Print("Z");
////    }
//    //  Render the scene
//    glFlush();
//    //  Make the rendered scene visible
//    glutSwapBuffers();
//}
//
///*
// *  GLUT calls this routine when an arrow key is pressed
// */
//void special(int key, int x, int y)
//{
//    //  Right arrow key - increase angle by 5 degrees
//    if (key == GLUT_KEY_RIGHT)
//        th += 5;
//    //  Left arrow key - decrease angle by 5 degrees
//    else if (key == GLUT_KEY_LEFT)
//        th -= 5;
//    //  Up arrow key - increase elevation by 5 degrees
//    else if (key == GLUT_KEY_UP)
//    {
//        if (ph + 5 < 90)
//        {
//            ph += 5;
//        }
//    }
//    //  Down arrow key - decrease elevation by 5 degrees
//    else if (key == GLUT_KEY_DOWN)
//    {
//        if (ph - 5 > 0)
//        {
//            ph -= 5;
//        }
//    }
//    //  Keep angles to +/-360 degrees
//    th %= 360;
//    ph %= 360;
//    //  Tell GLUT it is necessary to redisplay the scene
//    glutPostRedisplay();
//}
//
///*
// *  Set projection
// */
//void Project(double fov, double asp, double dim)
//{
//    //  Tell OpenGL we want to manipulate the projection matrix
//    glMatrixMode(GL_PROJECTION);
//    //  Undo previous transformations
//    glLoadIdentity();
//    //  Perspective transformation
//    if (fov)
//        gluPerspective(fov, asp, dim / 16, 16 * dim);
//    //  Orthogonal transformation
//    else
//        glOrtho(-asp * dim, asp * dim, -dim, +dim, -dim, +dim);
//    //  Switch to manipulating the model matrix
//    glMatrixMode(GL_MODELVIEW);
//    //  Undo previous transformations
//    glLoadIdentity();
//}
//
///*
// *  GLUT calls this routine when a key is pressed
// */
//void key(unsigned char ch, int x, int y)
//{
//    //  Exit on ESC
//    if (ch == 27)
//        exit(0);
//    //  Reset view angle
//    else if (ch == '0')
//        th = ph = 0;
//    //  Toggle axes
//    else if (ch == 'a' || ch == 'A')
//        axes = 1 - axes;
//    //  Change field of view angle
//    else if (ch == '-' && ch > 1)
//        fov++;
//    else if (ch == '=' && ch < 179)
//        fov--;
//    //  PageUp key - increase dim
//    else if (ch == GLUT_KEY_PAGE_DOWN) {
//        dim += 0.1;
//    }
//    //  PageDown key - decrease dim
//    else if (ch == GLUT_KEY_PAGE_UP && dim > 1) {
//        dim -= 0.1;
//    }
//    //  Keep angles to +/-360 degrees
//    th %= 360;
//    ph %= 360;
//    //  Reproject
//    Project(fov, asp, dim);
//    //  Tell GLUT it is necessary to redisplay the scene
//    glutPostRedisplay();
//}
//
///*
// *  GLUT calls this routine when the window is resized
// */
//void reshape(int width, int height)
//{
//    //  Ratio of the width to the height of the window
//    asp = (height > 0) ? (double)width / height : 1;
//    //  Set the viewport to the entire window
//    glViewport(0, 0, width, height);
//    //  Set projection
//    Project(fov, asp, dim);
//}
//
///*
// *  GLUT calls this toutine when there is nothing else to do
// */
//void idle()
//{
//    glutPostRedisplay();
//}
//
//int main(int argc, char* argv[])
//{
//    // Initialize GLUT and process user parameters
//    glutInit(&argc, argv);
//    // double buffered, true color 600*600
//    glutInitWindowSize(1000, 800);
//    glutInitDisplayMode(GLUT_RGB | GLUT_DEPTH | GLUT_DOUBLE);
//    // create the window
//    glutCreateWindow("Slight Spin_yj2563_cl3895");
//    //  Tell GLUT to call "idle" when there is nothing else to do
//    glutIdleFunc(idle);
//    //  Tell GLUT to call "display" when the scene should be drawn
//    glutDisplayFunc(display);
//    //  Tell GLUT to call "reshape" when the window is resized
//    glutReshapeFunc(reshape);
//    //  Tell GLUT to call "special" when an arrow key is pressed
//    glutSpecialFunc(special);
//    //  Tell GLUT to call "key" when a key is pressed
//    glutKeyboardFunc(key);
//    init();
//    //  Pass control to GLUT so it can interact with the user
//    glutMainLoop();
//    return 0;
//};
//
