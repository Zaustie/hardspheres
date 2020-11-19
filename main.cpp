#include <iostream>
#include <random>
#include <GL/freeglut.h>
#include <SDL2/SDL.h>
#include <cmath>
#include <iomanip>
// 320 Zeilen Darstellung

#define N 100 		// Number of Spheres
#define L 30.		// Size of the box
#define M 8 			// Number of subdivisions of the box
#define SIGMA 1.	// Radius of the spheres
#define DT 0.001	// Target time difference between two frames
#define T  5
#define MASS 1.

#define DT_MOV 0.01	// Target time difference between two frames
#define SCALE 0.1	// scale factor between simulation sizes and rendered size

#define EPSILON 1.00000001

#define OUTPUT_ON false

#define B1 10
#define B2 10

bool exitloop = false;
bool pause = true;

std::mt19937 rng;	// random number generator

// -----------------------------------------------------------------------------
// Quality of life functions

double rn_01()
{
	return double(rng())/double(rng.max());
}

double rn_gauss(double mu, double sigma)
{
	static double buffer = 0;
	static bool buffer_filled = false;

	double ret_value;

	if(buffer_filled)
	{
		ret_value = buffer;
		buffer_filled = false;
	}
	else
	{
		double u1 = rn_01();
		double u2 = rn_01();
		buffer = sqrt(-2*log(u1))*cos(2*M_PI*u2);
		ret_value = sqrt(-2*log(u1))*sin(2*M_PI*u2);
		buffer_filled = true;
	}
	return ret_value*sigma + mu;
}


void f_to_rgb(float in, float out[3])
{
	float color1[3] = {0.7f, 0.0f, 0.1f};
	float color2[3] = {0.7f, 0.7f, 0.7f};

	for(int k=0; k<3; ++k)
		out[k] = in*color1[k] + (1.-in)*color2[k];
}




// -----------------------------------------------------------------------------
// Simulation classes

class sphere_t;
class eventnode_t;

typedef sphere_t* cell_t;

class sphere_t
{
	public:
		double xp;
		double yp;
		double zp;
		double x;
		double y;
		double z;
		double vx;
		double vy;
		double vz;
		cell_t next;
};

enum event_enum
{
	event_none,
	event_cellcrossing,
	event_collision,
	event_timestep
};

class eventnode_t
{
	public:
		eventnode_t()
		{
			up = nullptr;
			left = nullptr;
			right = nullptr;
			idA = 0;
			idB = 0;
			time = 0;
			type = event_none;
		}
		eventnode_t* up;
		eventnode_t* left;
		eventnode_t* right;


		event_enum type;

		int idA;	// in 0 to N for cell crossing
		int idB;
		double time;
};

// -----------------------------------------------------------------------------
// Declerations of the following functions


eventnode_t* QuerryEvent(eventnode_t*);
void InsertEvent(eventnode_t*, eventnode_t*);
eventnode_t* GetEmptyEvent(eventnode_t*);
void RemoveEvent(eventnode_t*, eventnode_t*);
void AddEventToEmpty(eventnode_t*, eventnode_t*);

void CalculateCollisions(sphere_t*, cell_t, sphere_t*, eventnode_t*);
void RemoveCollisions(eventnode_t*, eventnode_t*, sphere_t*);
void ExecuteCollision(sphere_t*, sphere_t*);

// -----------------------------------------------------------------------------
// Init stuff

void InitSpheres(sphere_t* sphere, cell_t cell[M][M][M])
{
	int n=0;
	for(int i=0; i<M; ++i)
	{
		for(int j=0; j<M; ++j)
		{
			for(int k=0; k<M; ++k)
			{
				sphere[n].xp = L/M*(i+0.5);
				sphere[n].yp = L/M*(j+0.5);
				sphere[n].zp = L/M*(k+0.5);
				sphere[n].x = L/M*(i+0.5);
				sphere[n].y = L/M*(j+0.5);
				sphere[n].z = L/M*(k+0.5);
				cell[i][j][k] = sphere+n;
				++n;
				if(n==N)
					goto loop_end;
			}
		}
	}
	loop_end:

	if(N==8)
	{
		sphere[0].vx = 0.0000001;
		sphere[0].vy = 0.0000001;
		sphere[0].vz = 1.;
		sphere[0].next = nullptr;
		for(int i=1;i<N;++i)
		{
			sphere[i].vx = 0.0000001;
			sphere[i].vy = 0.0000001;
			sphere[i].vz = 0.0000001;
			sphere[i].next = nullptr;
		}

	}
	else
		for(int i=0; i<N; ++i)
		{
			sphere[i].vx = rn_gauss(0, sqrt(T));
			sphere[i].vy = rn_gauss(0, sqrt(T));
			sphere[i].vz = rn_gauss(0, sqrt(T));
			sphere[i].next = nullptr;
		}
}


void InitEventList(sphere_t* sphere, eventnode_t* eventnode, cell_t cell[M][M][M])
{
	eventnode_t* root = eventnode;
	eventnode_t* root_empty = eventnode + 1;
	root->time = 1e100;

	eventnode_t* frameevent = GetEmptyEvent(root_empty);
	frameevent->type = event_timestep;
	frameevent->time = DT;
	InsertEvent(frameevent, root);

	for(int i=0; i<N; ++i)
	{
		sphere_t* obj = sphere + i;

		int k = std::floor(obj->x/L*M);
		int l = std::floor(obj->y/L*M);
		int m = std::floor(obj->z/L*M);
		double tk = (obj->vx>0)?(((L/M)*(k+1)-obj->x)/obj->vx):(-(obj->x-k*(L/M))/obj->vx);
		double tl = (obj->vy>0)?(((L/M)*(l+1)-obj->y)/obj->vy):(-(obj->y-l*(L/M))/obj->vy);
		double tm = (obj->vz>0)?(((L/M)*(m+1)-obj->z)/obj->vz):(-(obj->z-m*(L/M))/obj->vz);
		double tmin = 0;

		if(tk<=tl && tk<=tm)
			tmin = tk;
		if(tl<=tk && tl<=tm)
			tmin = tl;
		if(tm<=tk && tm<=tl)
			tmin = tm;

		eventnode_t* newevent = GetEmptyEvent(root_empty);

		newevent->type = event_cellcrossing;
		newevent->time = tmin*EPSILON;
		newevent->idA = i;
		newevent->left = nullptr;
		newevent->right = nullptr;
		newevent->up = nullptr;

		InsertEvent(newevent, root);
	}

	for(int n=0; n<N; ++n)
	{
		for(int i=-1; i<2; ++i)
		{
			for(int j=-1; j<2; ++j)
			{
				for(int k=-1; k<2; ++k)
				{
					int a = std::floor(sphere[n].x/L*M);
					int b = std::floor(sphere[n].y/L*M);
					int c = std::floor(sphere[n].z/L*M);
					CalculateCollisions(sphere + n, cell[(a+i+M)%M][(b+j+M)%M][(c+k+M)%M], sphere, eventnode);
				}
			}
		}
	}
}

// -----------------------------------------------------------------------------
// Actual simulation stuff


eventnode_t* QuerryEvent(eventnode_t* root)
{
	if(root == nullptr)
	{
		std::cout << "QuerryEvent(eventnode_t*) needs a valid pointer!\n";
		exit(0);
	}


	eventnode_t* event = root;
	while(event->left!=nullptr)
	{
		event = event->left;
	}
	if(event->right != nullptr)
	{
		event->right->up = event->up;
		event->up->left = event->right;
	}
	else
	{
		event->up->left = nullptr;
	}
	event->up = nullptr;
	event->right = nullptr;
	event->left = nullptr;

	return event;
}


void InsertEvent(eventnode_t* event, eventnode_t* root)
{
	if(event == nullptr)
		return;

	eventnode_t* head = root;
	eventnode_t* next = nullptr;

	bool inserted = false;
	while(!inserted)
	{
		if(event->time < head->time)
		{
			next = head->left;
			if(next==nullptr)
			{
				head->left = event;
				event->up = head;
				inserted = true;
			}
		}
		else
		{
			next = head->right;
			if(next==nullptr)
			{
				head->right = event;
				event->up = head;
				inserted = true;
			}
		}
		head = next;
	}
}

eventnode_t* GetEmptyEvent(eventnode_t* root_empty)
{
	eventnode_t* ret = nullptr;
	if(root_empty->right != nullptr)
	{
		eventnode_t* tmp = root_empty->right;

		root_empty->right = tmp->right;
		root_empty->right->up = tmp->up;

		tmp->right = nullptr;
		tmp->up = nullptr;
		// make sure tmp is zero'd
		tmp->left = nullptr;
		tmp->idA  = 0;
		tmp->idB  = 0;
		tmp->time = 0;
		tmp->type = event_none;

		ret = tmp;
	}
	else
	{
		std::cout << "Fatal error in GetNewEvent(eventnode_t*)!" << std::endl;
		exit(0);
	}

	return ret;
}

void RemoveEvent(eventnode_t* event, eventnode_t* root)
{
	if(event->left != nullptr)
	{
		eventnode_t* tmp = event->left;
		tmp->up = event->up;
		if(event->up->left==event)
			event->up->left = tmp;
		else
			event->up->right = tmp;

		InsertEvent(event->right, root);
	}
	else if(event->right != nullptr)
	{
		eventnode_t* tmp = event->right;
		tmp->up = event->up;
		if(event->up->left==event)
			event->up->left = tmp;
		else
			event->up->right = tmp;
	}
	else
	{
		if(event->up->left==event)
			event->up->left = nullptr;
		else
			event->up->right = nullptr;
	}
}


void AddEventToEmpty(eventnode_t* event, eventnode_t* root_empty)
{
	eventnode_t* tmp = root_empty->right;
	tmp->up = event;
	event->right = tmp;
	event->up = root_empty;
	root_empty->right = event;
}



void HandleCellCrossing(eventnode_t* event, eventnode_t* eventnode, sphere_t* sphere, cell_t cell[M][M][M], double dt)
{
	sphere_t* obj = sphere + (event->idA);
	// delete obj from old cell
	int k1 = int(std::floor(obj->xp/L*M)+M)%M;
	int l1 = int(std::floor(obj->yp/L*M)+M)%M;
	int m1 = int(std::floor(obj->zp/L*M)+M)%M;
	cell_t tmp = cell[k1][l1][m1];
	if(cell[k1][l1][m1]==obj)
	{
		cell[k1][l1][m1] = obj->next;
		obj->next=nullptr;
	}
	else
	{
		while(tmp->next!=obj)
			tmp = tmp->next;
		tmp->next=obj->next;
		obj->next=nullptr;
	}

	// add obj to new cell
	int k2 = std::floor(obj->x/L*M);
	int l2 = std::floor(obj->y/L*M);
	int m2 = std::floor(obj->z/L*M);
	tmp = cell[k2][l2][m2];
	if(cell[k2][l2][m2]==nullptr)
	{
		cell[k2][l2][m2] = obj;
		obj->next = nullptr;
	}
	else
	{
		cell[k2][l2][m2] = obj;
		obj->next = tmp;
	}

	event->up = nullptr;
	event->right = nullptr;
	event->left = nullptr;

	// recalculate time for event and emplace into eventlist
	double tk = (obj->vx>0)?(((L/M)*(k2+1)-obj->x)/obj->vx):(-(obj->x-k2*(L/M))/obj->vx);
	double tl = (obj->vy>0)?(((L/M)*(l2+1)-obj->y)/obj->vy):(-(obj->y-l2*(L/M))/obj->vy);
	double tm = (obj->vz>0)?(((L/M)*(m2+1)-obj->z)/obj->vz):(-(obj->z-m2*(L/M))/obj->vz);
	double tmin = 0;
	if(tk<=tl && tk<=tm)
		tmin = tk;
	if(tl<=tk && tl<=tm)
		tmin = tl;
	if(tm<=tk && tm<=tl)
		tmin = tm;

	event->time = tmin*EPSILON;


// -----------------------------------------------------------------------------
// recalculating collisions
	RemoveCollisions(event, eventnode, sphere);

	// for all 9 new fields that were added, calculate new events.
	for(int i=-1; i<2; ++i)
	{
		for(int j=-1; j<2; ++j)
		{
			for(int k=-1; k<2; ++k)
			{
				CalculateCollisions(obj, cell[(k2+i+M)%M][(l2+j+M)%M][(m2+k+M)%M], sphere, eventnode);
			}
		}
	}

}


void CalculateCollisions(sphere_t* obj, cell_t cell_in, sphere_t* sphere, eventnode_t* eventnode)
{
	eventnode_t* root = eventnode;
	eventnode_t* root_empty = eventnode + 1;
	sphere_t* current = cell_in;
	while(current != nullptr)
	{
		int idA = std::min(current-sphere, obj-sphere);
		int idB = std::max(current-sphere, obj-sphere);

		bool isinlist = false;
		for(int i=2; i<N*(N+1)+2+2; ++i)
		{
			if(eventnode[i].type == event_collision && idA==eventnode[i].idA && idB==eventnode[i].idB)
			{
				isinlist = true;
				break;
			}
		}

		if(!isinlist)
		{
			double s1[3] = {obj->x,  obj->y,  obj->z };
			double v1[3] = {obj->vx, obj->vy, obj->vz};
			double s2[3] = {current->x,  current->y,  current->z };
			double v2[3] = {current->vx, current->vy, current->vz};
			double Ds[3] = {0};

			for(int i=0; i<3; ++i)
			{
				if(s1[i]<s2[i])
				{
					if(fabs(s2[i]-s1[i]) > fabs(s2[i]-L-s1[i]))
						if(s1[i]<L/M or s2[i]>L-L/M)
							s2[i] -= L;
				}
				else
				{
					if(fabs(s2[i]-s1[i]) > fabs(s2[i]-(s1[i]-L)))
						if(s2[i]<L/M or s1[i]>L-L/M)
							s1[i] -= L;
				}
				Ds[i] = s2[i]-s1[i];
			}


			double a = 0;
			double b = 0;
			double c = 0;

			for(int i=0; i<3; ++i)
			{
				a += (v2[i]-v1[i])*(v2[i]-v1[i]);
				b += (v2[i]-v1[i])*Ds[i];
				c += Ds[i]*Ds[i];
			}
			c -= 4*SIGMA*SIGMA;

			if(b*b-a*c>0)
			{
				//Collision happens and we need to add a new event to the tree
				double t1 = (-b-sqrt(b*b-a*c))/a;
				double t2 = (-b+sqrt(b*b-a*c))/a;
				double t = -1;
				if(t1>1e-14)
					t = t1;
//				else if(t2>1e-14)
// 					t = t2;
				if(t>1e-14)
				{
					eventnode_t* newevent = GetEmptyEvent(root_empty);

					newevent->time = t/EPSILON; // /EPSILON;
					newevent->idA = std::min(int(obj-sphere), int(current-sphere));
					newevent->idB = std::max(int(obj-sphere), int(current-sphere));
					newevent->type = event_collision;

					InsertEvent(newevent, root);
				}
			}

		}
		current = current->next;
	}
}



void RemoveCollisions(eventnode_t* event, eventnode_t* eventnode, sphere_t* sphere)
{
	eventnode_t* root = eventnode;
	eventnode_t* root_empty = eventnode+1;

	for(int i=3; i<N*(N+1)+2+2; ++i)
	{
		if((root[i].idA == event->idA or root[i].idB == event->idA) && root[i].type == event_collision)
		{
			RemoveEvent(root+i, root);
			root[i].type = event_none;
			root[i].idA = 0;
			root[i].idB = 0;
			root[i].time = 0;
			root[i].left = nullptr;
			root[i].up = nullptr;
			root[i].right = nullptr;
			AddEventToEmpty(root+i, root_empty);
		}
	}
}

void ExecuteCollision(sphere_t* s1, sphere_t* s2)
{
	double pos1[3] = {s1->x,  s1->y,  s1->z};
	double pos2[3] = {s2->x,  s2->y,  s2->z};
	double vel1[3] = {s1->vx, s1->vy, s1->vz};
	double vel2[3] = {s2->vx, s2->vy, s2->vz};
	double m1 = MASS;
	double m2 = MASS;


	for(int i=0; i<3; ++i)
	{

		if(pos1[i]<pos2[i])
		{
			if(fabs(pos2[i]-pos1[i]) > fabs(pos2[i]-L-pos1[i]))
				if(pos1[i]<L/M or pos2[i]>L-L/M)
					pos2[i] -= L;
		}
		else
		{
			if(fabs(pos2[i]-pos1[i]) > fabs(pos2[i]-(pos1[i]-L)))
				if(pos2[i]<L/M or pos1[i]>L-L/M)
					pos1[i] -= L;
		}
	}

	double k[3] = {0};
	double a = 0;

	double dist = 0;
	for(int i=0; i<3; ++i)
	{
		k[i] = pos2[i] - pos1[i];
		dist += k[i]*k[i];
	}
	dist = sqrt(dist);
	for(int i=0; i<3; ++i)
		k[i] /= dist;

	for(int i=0; i<3; ++i)
		a += k[i]*(vel2[i]-vel1[i]);
	a *= 2./(1./m1 + 1./m2);


	s1->vx += a/m1*k[0];
	s1->vy += a/m1*k[1];
	s1->vz += a/m1*k[2];

	s2->vx -= a/m2*k[0];
	s2->vy -= a/m2*k[1];
	s2->vz -= a/m2*k[2];
}



// -----------------------------------------------------------------------------
// -------------------------Graphics stuff starts here--------------------------
// -----------------------------------------------------------------------------

namespace cam
{
	const float rotate_a_set = 160.0f;
	const float acc_set = 3.2f*SCALE*L;
	float rotate[3] = {0.0f, 0.0f, 0.0f};
	float rotate_v[3] = {0.0f, 0.0f, 0.0f};
	int mom_active[3] = {0, 0, 0};
	float pos[3] = {0.0f, 0.0f, -SCALE*L*1.4};
	float vel[3] = {0.0f, 0.0f, 0.0f};
	int acc_active[3] = {0, 0, 0};
	float dt = DT_MOV;
}

sphere_t* sphere_draw_ptr = nullptr;

void Draw_Spheres ()
{
    glColor3f (1.0, 1.0, 1.0);
	glPushMatrix();
		glutWireCube(SCALE*L);
	glPopMatrix();


	if(sphere_draw_ptr!=nullptr)
	{
		float distance[N];
		for(int i=0; i<N; ++i)
			distance[i] = (SCALE*(sphere_draw_ptr[i].x-0.5*L)+cam::pos[0])
							*(SCALE*(sphere_draw_ptr[i].x-0.5*L)+cam::pos[0])
						+ (SCALE*(sphere_draw_ptr[i].y-0.5*L)+cam::pos[1])
							*(SCALE*(sphere_draw_ptr[i].y-0.5*L)+cam::pos[1])
						+ (SCALE*(sphere_draw_ptr[i].z-0.5*L)+cam::pos[2])
							*(SCALE*(sphere_draw_ptr[i].z-0.5*L)+cam::pos[2]);
		int order[N];
		for(int i=0; i<N; ++i)
			order[i] = i;

		for(int i=0; i<N; ++i)
		{
			float distance_max = distance[i];
			int index = i;
			for(int k=i+1; k<N; ++k)
			{
				if(distance_max<distance[k])
				{
					distance_max = distance[k];
					index = k;
				}
			}
			float tmp_distance = distance[i];
			distance[i] = distance_max;
			distance[index] = tmp_distance;
			int tmp_index = order[i];
			order[i] = order[index];
			order[index] = tmp_index;
		}

		float minspeed = 100000;
		float maxspeed = 0;
		float speed = 0;
		for(int i=0; i<N; ++i)
		{
			speed = sphere_draw_ptr[i].vx*sphere_draw_ptr[i].vx
				  + sphere_draw_ptr[i].vy*sphere_draw_ptr[i].vy
				  + sphere_draw_ptr[i].vz*sphere_draw_ptr[i].vz;
			if(speed > maxspeed)
				maxspeed = speed;
			if(minspeed > speed)
				minspeed = speed;
		}

		for(int i=0; i < N; i++)  
		{
			float color[3] = {0.9, 0.6, 0.6};

			speed = sphere_draw_ptr[order[i]].vx*sphere_draw_ptr[order[i]].vx
				  + sphere_draw_ptr[order[i]].vy*sphere_draw_ptr[order[i]].vy
				  + sphere_draw_ptr[order[i]].vz*sphere_draw_ptr[order[i]].vz;
			f_to_rgb((speed-minspeed)/(maxspeed-minspeed), color);

			glColor3f(color[0], color[1], color[2]);
			glPushMatrix();
				glTranslatef (SCALE*((float)sphere_draw_ptr[order[i]].x - 0.5*L),
							  SCALE*((float)sphere_draw_ptr[order[i]].y - 0.5*L),
							  SCALE*((float)sphere_draw_ptr[order[i]].z - 0.5*L));
				glShadeModel(GL_SMOOTH);
//				glShadeModel(GL_FLAT);
				glutSolidSphere(SCALE*SIGMA, B1, B2);
//				glutWireSphere(SIGMA, B1, B2);
			glPopMatrix();
		}
	}
}

void CalcCamPos(void)
{
	GLfloat matrix[16];
	glGetFloatv(GL_MODELVIEW_MATRIX, matrix);
	float dreibein[3][3];
	for(int i=0; i<3; ++i)
		for(int k=0; k<3; ++k)
			dreibein[k][i] = matrix[4*i + k];
	// dreibein[0] gibt den von der Kamera nach rechts gehenden Vektor (dos)
	// dreibein[1] gibt den von der Kamera nach oben gehenden Vektor (dov)
	// dreibein[2] gibt den von der Kamera nach vorne gehenden Vektor (dof)

	float acc[3] = {0.0f, 0.0f, 0.0f};
	for (int k=0; k<3; ++k)
		for (int i=0; i<3; ++i)
			acc[i] += cam::acc_set*cam::acc_active[k]*dreibein[k][i];

	for(int i=0; i<3; ++i)
	{
		cam::vel[i] += cam::dt*acc[i];	
		cam::pos[i] += cam::dt*cam::vel[i];
		cam::vel[i]  = std::min(cam::vel[i], cam::acc_set);
		cam::vel[i]  = std::max(cam::vel[i], -cam::acc_set);
	}
	for(int i=0; i<3; ++i)
	{
		cam::vel[i] *= 0.9;
		if(cam::vel[i]<0.001 && cam::vel[i]>0.001)
			cam::vel[i] = 0.0;
	}


	float rotate_a[3] = {0.0f, 0.0f, 0.0f};
	for (int i=0; i<3; ++i)
		rotate_a[i] -= cam::rotate_a_set*cam::mom_active[i];

	for(int i=0; i<3; ++i)
	{
		cam::rotate_v[i] += cam::dt*rotate_a[i];	
		cam::rotate[i] += cam::dt*cam::rotate_v[i];
		cam::rotate_v[i]  = std::min(cam::rotate_v[i], cam::rotate_a_set);
		cam::rotate_v[i]  = std::max(cam::rotate_v[i], -cam::rotate_a_set);
	}
	for(int i=0; i<3; ++i)
	{
		cam::rotate_v[i] *= 0.9;
		if(cam::rotate_v[i]<0.001 && cam::rotate_v[i]>0.001)
			cam::rotate_v[i] = 0.0;
	}
}


void Init()
{
	GLfloat sun_direction[] = { 10.0f, 10.0f, -1.0f, 0.0f};
	GLfloat sun_intensity[] = { 0.7f, 0.7f, 0.7f, 1.0f};
	GLfloat ambient_intensity[] = { 0.3f, 0.3f, 0.3f, 1.0f};

	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

    glMatrixMode(GL_MODELVIEW);

	glEnable(GL_DEPTH_TEST);            // Draw only closest surfaces

	glEnable(GL_NORMALIZE);
	glEnable(GL_CULL_FACE);


	glEnable(GL_LIGHTING);              // Set up ambient light.
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient_intensity);

	glEnable(GL_LIGHT0);                // Set up sunlight.
	glLightfv(GL_LIGHT0, GL_POSITION, sun_direction);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, sun_intensity);

	glEnable(GL_COLOR_MATERIAL);        // Configure glColor().
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);

}

void display (void)
{
	GLfloat sun_direction[] = { 10.0f, 10.0f, -1.0f, 0.0f};

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glLoadIdentity(); // Load the Identity Matrix to reset drawing locations


	glRotatef(cam::rotate[0], 1.0,0.0,0.0);
	glRotatef(cam::rotate[1], 0.0,1.0,0.0);
	glRotatef(cam::rotate[2], 0.0,0.0,1.0);
	glTranslatef(cam::pos[0], cam::pos[1], cam::pos[2]);

	glLightfv(GL_LIGHT0, GL_POSITION, sun_direction);

	CalcCamPos();

	Draw_Spheres();

	glutSwapBuffers();
}



void reshape(int width, int height)
{
	glViewport(0, 0, (GLsizei)width, (GLsizei)height);
	glMatrixMode(GL_PROJECTION);	// Switch to the projection matrix so that we can manipulate how our scene is viewed
	glLoadIdentity();				// Reset the projection matrix to the identity matrix so that we don't get any artifacts (cleaning up) 
	gluPerspective(60, (GLfloat)width / (GLfloat)height, 0.01, 500.0); // Set the Field of view angle (in degrees), the aspect ratio of our window, and the new and far planes
	glMatrixMode(GL_MODELVIEW); // Switch back to the model view matrix, so that we can start drawing shapes correctly
}

static void KeyUp(unsigned char key, int x, int y)
{
    switch (key) {
		case 'w':
			cam::acc_active[2] = std::max(cam::acc_active[2]-1, -1);
			glutPostRedisplay();
			break;
		case 's':
			cam::acc_active[2] = std::min(cam::acc_active[2]+1, 1);
			glutPostRedisplay();
			break;
		case 'a':
			cam::acc_active[0] = std::max(cam::acc_active[0]-1, -1);
			glutPostRedisplay();
			break;
		case 'd':
			cam::acc_active[0] = std::min(cam::acc_active[0]+1, 1);
			glutPostRedisplay();
			break;
		case 'e':
			cam::acc_active[1] = std::min(cam::acc_active[1]+1, 1);
			glutPostRedisplay();
			break;
		case 'q':
			cam::acc_active[1] = std::max(cam::acc_active[1]-1, -1);
			glutPostRedisplay();
			break;
      case 27:
		exit(0);
    }
}


static void Key(unsigned char key, int x, int y)
{
    switch (key) {
		case 'w':
			cam::acc_active[2] = std::min(cam::acc_active[2]+1, 1);
			glutPostRedisplay();
			break;
		case 's':
			cam::acc_active[2] = std::max(cam::acc_active[2]-1, -1);
			glutPostRedisplay();
			break;
		case 'a':
			cam::acc_active[0] = std::min(cam::acc_active[0]+1, 1);
			glutPostRedisplay();
			break;
		case 'd':
			cam::acc_active[0] = std::max(cam::acc_active[0]-1, -1);
			glutPostRedisplay();
			break;
		case 'e':
			cam::acc_active[1] = std::max(cam::acc_active[1]-1, -1);
			glutPostRedisplay();
			break;
		case 'q':
			cam::acc_active[1] = std::min(cam::acc_active[1]+1, 1);
			glutPostRedisplay();
			break;
		case ' ':
			pause ^= true;
			glutPostRedisplay();
			break;
      case 27:	// "esc"
			exitloop = true;
			//exit(0);
			break;
    }
}

static void SpecialKeyUp(int key, int x, int y)
{
    switch (key) {
		case GLUT_KEY_UP:
			glutPostRedisplay();
			cam::mom_active[0] = std::max(cam::mom_active[0]-1, -1);
			break;
		case GLUT_KEY_DOWN:
			cam::mom_active[0] = std::min(cam::mom_active[0]+1, 1);
			glutPostRedisplay();
			break;
		case GLUT_KEY_LEFT:
			cam::mom_active[1] = std::max(cam::mom_active[1]-1, -1);
			glutPostRedisplay();
			break;
		case GLUT_KEY_RIGHT:
			cam::mom_active[1] = std::min(cam::mom_active[1]+1, 1);
			glutPostRedisplay();
			break;
    }
}

static void SpecialKey(int key, int x, int y)
{
    switch (key) {
		case GLUT_KEY_UP:
			cam::mom_active[0] = std::min(cam::mom_active[0]+1, 1);
			glutPostRedisplay();
			break;
		case GLUT_KEY_DOWN:
			cam::mom_active[0] = std::max(cam::mom_active[0]-1, -1);
			glutPostRedisplay();
			break;
		case GLUT_KEY_LEFT:
			cam::mom_active[1] = std::min(cam::mom_active[1]+1, 1);
			glutPostRedisplay();
			break;
		case GLUT_KEY_RIGHT:
			cam::mom_active[1] = std::max(cam::mom_active[1]-1, -1);
			glutPostRedisplay();
			break;
    }
}


// --------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{

// -----------------------------------------------------------------------------
// Graphics preparation
	glutInit(&argc, argv);				// Initialize GLUT
	glutInitDisplayMode(GLUT_SINGLE);	// Set up a basic display buffer (only single buffered for now)
	glutInitWindowSize(2000, 2000);		// Set the width and height of the window
	glutInitWindowPosition(800, 200);	// Set the position of the window
	glutCreateWindow("Spheres");		// Set the title for the window

	Init(); 						// Init some stuff for rendering
	glutDisplayFunc (display);		// rendering function
	glutIdleFunc    (display);		// idle funciton
	glutReshapeFunc (reshape);		// window reshaping function
	glutKeyboardFunc(Key);			// regular keypress function (ascii letters)
	glutKeyboardUpFunc(KeyUp);
	glutSpecialFunc(SpecialKey);	// special keypress function (arrows, etc.…)
	glutSpecialUpFunc(SpecialKeyUp);

	glutMainLoopEvent();			// Run one iteration of gluts MainLoop


// -----------------------------------------------------------------------------
// Setup for simulation + main loop 
	if(OUTPUT_ON)
		std::cout << std::setprecision(8);

	eventnode_t* eventnode = new eventnode_t[N*(N+1)+2+2];	// 0 is root, 1 is root of emptys
	sphere_t* sphere = new sphere_t[N];
	sphere_draw_ptr = sphere;
	// maximum N*N collisions, N cellcrossings, 1 frameevent, 1 end event, 1 root, 1 root of emptys.

	for(int i=1; i<N*(N+1)+1; ++i)
	{
		eventnode[i].right = eventnode + i+1;
		eventnode[i+1].up = eventnode + i;
	}

	cell_t cell[M][M][M] = {0};

	InitSpheres(sphere, cell);

	InitEventList(sphere, eventnode, cell);

	double t = 0.;
	double dt = DT;
	while(!exitloop)
	{
		if(!pause)
		{
			// reduce time of all events in eventtree by next event's time
			eventnode_t* event = QuerryEvent(eventnode);
			dt = event->time;

			eventnode_t* tmp = event;
			for(int i=2; i<N*(N+1)+2+2; ++i)
			{
				if(eventnode[i].type != event_none)
				{
					eventnode[i].time -= dt;
				}
			}

			for(int i=0; i<N; ++i)
			{
				sphere[i].xp = sphere[i].x;
				sphere[i].yp = sphere[i].y;
				sphere[i].zp = sphere[i].z;
				sphere[i].x += dt*sphere[i].vx;
				sphere[i].y += dt*sphere[i].vy;
				sphere[i].z += dt*sphere[i].vz;
			}
			for(int i=0; i<N; ++i)
			{
				if(sphere[i].x>=L)
					sphere[i].x -= L;
				if(sphere[i].x<0)
					sphere[i].x += L;
				if(sphere[i].y>=L)
					sphere[i].y -= L;
				if(sphere[i].y<0)
					sphere[i].y += L;
				if(sphere[i].z>=L)
					sphere[i].z -= L;
				if(sphere[i].z<0)
					sphere[i].z += L;
			}
			switch(event->type)
			{
				case event_cellcrossing:
				{
//std::cout << "cellcrossing\n";
					HandleCellCrossing(event, eventnode, sphere, cell, dt);

					InsertEvent(event, eventnode);

					break;
				}
				case event_collision:
				{
//std::cout << "collision\n";
					// execute collision and recalculate new events for both participating particles
					int idA = event->idA;
					int idB = event->idB;

					ExecuteCollision(sphere + idA, sphere + idB);


					for(int i=3; i<3+N; ++i)
					{
						if(eventnode[i].type == event_cellcrossing && (event->idA == eventnode[i].idA or event->idB == eventnode[i].idA))
						{
							//Recalc cell crossing events
							sphere_t* obj = sphere + eventnode[i].idA;

							int k2 = std::floor(obj->x/L*M);
							int l2 = std::floor(obj->y/L*M);
							int m2 = std::floor(obj->z/L*M);
							double tk = (obj->vx>0)?(((L/M)*(k2+1)-obj->x)/obj->vx):(-(obj->x-k2*(L/M))/obj->vx);
							double tl = (obj->vy>0)?(((L/M)*(l2+1)-obj->y)/obj->vy):(-(obj->y-l2*(L/M))/obj->vy);
							double tm = (obj->vz>0)?(((L/M)*(m2+1)-obj->z)/obj->vz):(-(obj->z-m2*(L/M))/obj->vz);
							double tmin = 0;
							if(tk<=tl && tk<=tm)
								tmin = tk;
							if(tl<=tk && tl<=tm)
								tmin = tl;
							if(tm<=tk && tm<=tl)
								tmin = tm;


							RemoveEvent(eventnode + i, eventnode);
							eventnode[i].time = tmin*EPSILON;
							eventnode[i].up = nullptr;
							eventnode[i].left = nullptr;
							eventnode[i].right = nullptr;
							InsertEvent(eventnode + i, eventnode);
						}
					}

					event->type = event_none;
					event->idA = 0;
					event->idB = 0;
					event->time = 0;
					AddEventToEmpty(event, eventnode + 1);

					for(int i=2; i<N*(N+1)+2+2; ++i)
					{
						if(eventnode[i].type == event_collision && (event->idA == eventnode[i].idA or event->idA == eventnode[i].idB or event->idB == eventnode[i].idA or event->idB == eventnode[i].idB))
						{
							RemoveEvent(eventnode + i, eventnode);
							eventnode[i].type = event_none;
							eventnode[i].idA = 0;
							eventnode[i].idB = 0;
							eventnode[i].time = 0;
							AddEventToEmpty(eventnode + i, eventnode + 1);
						}
					}

					for(int i=-1; i<2; ++i)
					{
						for(int j=-1; j<2; ++j)
						{
							for(int k=-1; k<2; ++k)
							{
								int a = std::floor(sphere[idA].x/L*M);
								int b = std::floor(sphere[idA].y/L*M);
								int c = std::floor(sphere[idA].z/L*M);
								CalculateCollisions(sphere + idA, cell[(a+i+M)%M][(b+j+M)%M][(c+k+M)%M], sphere, eventnode);
							}
						}
					}
					for(int i=-1; i<2; ++i)
					{
						for(int j=-1; j<2; ++j)
						{
							for(int k=-1; k<2; ++k)
							{
								int a = std::floor(sphere[idB].x/L*M);
								int b = std::floor(sphere[idB].y/L*M);
								int c = std::floor(sphere[idB].z/L*M);
								CalculateCollisions(sphere + idB, cell[(a+i+M)%M][(b+j+M)%M][(c+k+M)%M], sphere, eventnode);
							}
						}
					}

					break;
				}
				case event_timestep:
				{
//std::cout << "timestep\n";
					event->time = DT;
					InsertEvent(event, eventnode);
					if(OUTPUT_ON)
					{
						t += DT;
						std::cout << t;
						for(int i=0; i<N; ++i)
							std::cout << " " << sqrt(sphere[i].vx*sphere[i].vx + 
													 sphere[i].vy*sphere[i].vy +
												 	 sphere[i].vz*sphere[i].vz);
						std::cout << "\n";
					}
					glutPostRedisplay();
					glutMainLoopEvent();
					// TODO: return frametime
					break;
				}
				default:
					glutPostRedisplay();
					glutMainLoopEvent();
					std::cout << "NANI?\n";
					exitloop = true;
					break;
			}
		}
		else
		{
			glutPostRedisplay();
			glutMainLoopEvent();
			// TODO: return frametime
		}
	}

	delete[] sphere;
	delete[] eventnode;

	return 0;
}
