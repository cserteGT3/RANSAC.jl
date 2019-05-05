// Declarations


// Code

/* filter variables*/
static float z_1=0.0, z_2=0.0, z_3=0.0;
static float ztmp_1=0.0, ztmp_2=0.0;

/* filter parameters: LPF3*/
float ad11 = 0.99590939195495;
float ad12 = 0.0009598654638844413;
float ad13 = 3.652013455242641e-7;
float ad21 = -11.32353396334311;
float ad22 = 0.8877776077432195;
float ad23 = 0.0006156713046306184;
float ad31 = -19089.67481549373;
float ad32 = -193.61650049754638;
float ad33 = 0.30752107344727025;
float bd1 = 0.004090608044583261;
float bd2 = 11.323533963343152;
float bd3 = 19089.67481549386;

ztmp_1 = ad11*z_1 + ad12*z_2 + ad13*z_3 + bd1*data.velocity;
ztmp_2 = ad21*z_1 + ad22*z_2 + ad23*z_3 + bd2*data.velocity;
z_3 = ad31*z_1 + ad32*z_2 + ad33*z_3 + bd3*data.velocity;

z_1 = ztmp_1;
z_2 = ztmp_2;

data.filt_vel = z_1;

int start = 0;
int ampl = 30;
int period = 1500;
if (ticks < start)
{data.torque = 0;}
else
{data.torque = (((ticks-start)%period-period/2) > 0 ? ampl : 0);}