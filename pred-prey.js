// Lotke-Volterra predator-prey system ODE:
//   dx/dt = ax - bxy = f(x,y)
//   dy/dt = dxy - cy = g(x,y)
//   @Params: a,b,c,d,xo,yo > 0, all reals
function Model (a,b,c,d,xo,yo) {
	console.log("Model: a:" + a + " b:" + b + " c:" + c + " d:" + d + " xo:" + xo + " yo:" + yo);
	this.t = 0;
	this.a = a;
	this.b = b;
	this.c = c;
	this.d = d;
	this.x = xo;
	this.y = yo;
	this.f = function(x,y) {
		return this.a * x - this.b * x * y;
	};
	this.g = function(x, y) {
		return this.d * x * y - this.c * y;
	};
	// rk4 step function
	// @param h is timestep
	this.step = function(h) {
		// console.log("Model: step: x:" + this.x + " y:" + this.y);
		var k1 = h * this.f(this.x, this.y);
		var l1 = h * this.g(this.x, this.y);

		var k2 = h * this.f(this.x + k1 / 2, this.y + l1 / 2);
		var l2 = h * this.g(this.x + k1 / 2, this.y + l1 / 2);

		var k3 = h * this.f(this.x + k2 / 2, this.y + l2 / 2);
		var l3 = h * this.g(this.x + k2 / 2, this.y + l2 / 2);

		var k4 = h * this.f(this.x + k3, this.y + l3);
		var l4 = h * this.g(this.x + k3, this.y + l3);

		var k = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
		var l = (l1 + 2 * l2 + 2 * l3 + l4) / 6;

		this.t += h;
		this.x += k;
		this.y += l;
		if (this.x > 1) this.x = 1;
		if (this.y > 1) this.y = 1;
	};
};