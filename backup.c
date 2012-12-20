int simple_orbit_runge_kutta_4th(double * vars_in, double * vars_out, double step){
	/*
	Our variables are as follows:
	0 --> x-position of satellite
	1 --> y-position of satellite
	2 --> x-velocity of satellite
	3 --> y-velocity of satellite
	4 --> time
	CONSTANTS:
	5 --> mass of object
	6 --> time limit
	USEFUL THINGS TO PRINT:
	7 --> total energy / unit mass of satellite
	*/

	// copy the constants first
	int i;
	for(i=5;i<8;i++){
		vars_out[i] = vars_in[i];
	}

	// there are 4 variables, and this is a 4th order method, so we need 4 * 4 = 16
	// k values.

	//  			   |0--xpos---|  |1--ypos---|  |2--xvel---|  |3--yvel---|
	double k[4][4] = { {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0} };

	// this is slightly confusing because my indices are the other way 
	// round to how they are written in the lecture notes. I like my indices
	// better.

	double vel_derivative_multiplier = 
			    -vars_in[5]*GRAVITATIONAL_CONSTANT/
				pow(sqrt(pow(vars_in[0], 2) + pow(vars_in[1], 2)),3);

	// k1
	k[0][0] = vars_in[2];
	k[1][0] = vars_in[3];

	k[2][0] = - vars_in[0]*vel_derivative_multiplier;
	k[3][0] = - vars_in[1]*vel_derivative_multiplier;

	// k2
	k[0][1] = vars_in[2] + step * k[2][0]/2;
	k[1][1] = vars_in[3] + step * k[3][0]/2;

	k[2][1] = vel_derivative_multiplier * pow(vars_in[2] + step * k[2][0]/2, 2);
	k[3][1] = vel_derivative_multiplier * pow(vars_in[3] + step * k[3][0]/2, 2);

	// k3
	k[0][2] = vars_in[2] + step * k[2][1] / 2;
	k[1][2] = vars_in[3] + step * k[3][1] / 2;

	k[2][2] = vel_derivative_multiplier * pow(vars_in[2] + step * k[2][1] / 2, 2);
	k[3][2] = vel_derivative_multiplier * pow(vars_in[3] + step * k[3][1] / 2, 2);

	// k4
	k[0][3] = vars_in[2] * step * k[2][2];
	k[1][3] = vars_in[3] * step * k[3][2];

	k[2][3] = vel_derivative_multiplier * (vars_in[2] + step * k[2][2]);
	k[3][3] = vel_derivative_multiplier * (vars_in[3] + step * k[3][2]);


	for(i=0;i<4;i++){
		vars_out[i] = vars_in[i] + step/6*(k[i][0] + 2*k[i][1] + 2*k[i][2] + k[i][3]);
	}

	vars_out[4] = vars_in[4] + step;
	vars_out[7] = .5 * (pow(vars_out[2], 2) - pow(vars_out[3], 2)) + 
		GRAVITATIONAL_CONSTANT * vars_out[5] / sqrt(pow(vars_out[0],2) + pow(vars_out[1],2));

	// check we have not exceeded the time limit:
	if(vars_out[4] >= vars_in[6]){
		return STOP_ITERATING;
	}else{
		return CONTINUE_ITERATING;
	}

}