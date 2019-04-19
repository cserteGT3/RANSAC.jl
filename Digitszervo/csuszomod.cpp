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

data.pos_err = pos_ref + init_pos - data.position;
data.vel_err = -1*data.velocity;
data.sigma = data.pos_err + lambda*data.vel_err;

if(data.sigma < 0){
	data.torque = -1*M;}
else if(data.sigma > 0){
	data.torque = M;}
else{data.torque = 0;}
