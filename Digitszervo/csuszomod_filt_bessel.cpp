// Declarations

sigma
pos_err
vel_err

// Code

float pos_ref = 99.903;
float lambda = 2;
float M = 20;
float delta_10rad = 814.8733;

static double init_pos = 0;
if (ticks==1){init_pos = data.position + delta_10rad;}


/* filter variables*/
static float z_1=0.0, z_2=0.0, z_3=0.0;
static float ztmp_1=0.0, ztmp_2=0.0;

/* filter parameters: Bessel3*/
float ad11 = 0.9519298974843171;
float ad12 = 0.0008337130255001748;
float ad13 = 2.600911144761245e-7;
float ad21 = -120.96685586301366;
float ad22 = 0.5668804362636986;
float ad23 = 0.00034345282479987223;
float ad31 = -159737.89968605863;
float ad32 = -629.4283825447227;
float ad33 = -0.08051288648132714;
float bd1 = 0.04807010251378863;
float bd2 = 120.96685586292867;
float bd3 = 159737.89968584484;

ztmp_1 = ad11*z_1 + ad12*z_2 + ad13*z_3 + bd1*data.velocity;
ztmp_2 = ad21*z_1 + ad22*z_2 + ad23*z_3 + bd2*data.velocity;
z_3 = ad31*z_1 + ad32*z_2 + ad33*z_3 + bd3*data.velocity;

z_1 = ztmp_1;
z_2 = ztmp_2;

data.velocity = z_1;

data.pos_err = pos_ref + init_pos - data.position;
data.vel_err = -1*data.velocity;
data.sigma = data.pos_err + lambda*data.vel_err;

if(data.sigma < 0){
	data.torque = -1*M;}
else if(data.sigma > 0){
	data.torque = M;}
else{data.torque = 0;}
