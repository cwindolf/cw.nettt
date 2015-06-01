// from http://en.wikipedia.org/wiki/Shallow_water_equations
// conservative form of shallow water equations:
// for 		n		- 	fluid column height at x,y,t
//			(u,v)	-	fluid horizontal velocity at x,y,t
//			g		-	acceleration due to gravity
//			H 		- 	mean water height
// let's write the equations down:
// if 		   	{  n  }
//		U  := 	{ n*u }
//			   	{ n*v }
//
//				{ n*u               }
//		F(U) :=	{ n*u^2 + g*n^2 / 2 }
//				{ n*u*v             }
//
//				{ n*v }
//		G(U) := { n*u*v}
//				{ n*v^2 + g*n^2 / 2 }
// then we can write the system as:
//
//		dU     dF(U)     dG(U)     _
//		--  +  -----  +  -----  =  0
//		dt      dx        dy
//
// use Lax-Wendroff scheme per https://www.mathworks.com/moler/exm/chapters/water.pdf
// at any given moment, we have 2d array U_n, which has values at integer points (x,y)
// we want U_n+1. won't write out the method, but it's implemented below.

/* ********************************************************************************* */
// class ShallowWater: interactive numerical shallow water system
// takes parameters above.
// reflective BCs: u = 0 at y = 0,QUANT
//				   v = 0 at x = 0,QUANT
// so our parameters are then:
//		QUANT		-	we're doing numerics. so we need a number of discrete points
//				    	to sample. our following parameters (except H) will thus be
//						QUANTxQUANT arrays.
//		LENGTH		-   the physical side length of the bathtub
// 		n_o			-	QUANTxQUANT float array. initial value for column height.
//		u_o, v_o	-	same. IV for velocity. must conform to reflective BCs above.
function SWE(QUANT, LENGTH, n_o, u_o, v_o, g) {
	this.QUANT = QUANT || 100;
	this.LENGTH = LENGTH || 10;
	this.dd = this.LENGTH / this.QUANT;
	// if n_o is not supplied, let's make a barely filled bathtub.
	this.n;
	if (n_o == undefined) {
		console.log("create n");
		this.n = [];
		for (var _i = 0; _i < QUANT; _i++) {
			this.n.push([]);
			for (var _j = 0; _j < QUANT; _j++) {
				this.n[_i].push(1);
			}
		}
	} else {
		this.n = n_o;
	}
	this.u = u_o || this.n.map(function(a) { return a.map(function() { return 0; }) }); // backup is a correct-sized
	this.v = v_o || this.n.map(function(a) { return a.map(function() { return 0; }) }); // array of zeros
	this.g = g || -9.807; // meters per second^2
	this.H = 1; // calculate the average water height in the IC

	// System
	this.U = function(x, y) {
		if (x < 0) x = 0;
		if (y < 0) y = 0;
		if (x >= this.QUANT) x = this.QUANT - 1;
		if (y >= this.QUANT) y = this.QUANT - 1;
		return [this.n[x][y], 
				this.n[x][y] * this.u[x][y], 
				this.n[x][y] * this.v[x][y]];
	};
	this.F = function(U_xy) {
		return [U_xy[0] * U_xy[1], 
				U_xy[0] * U_xy[1] * U_xy[1] + this.g * U_xy[0] * U_xy[0] / 2,
				U_xy[0] * U_xy[1] * U_xy[2]];
	};
	this.G = function(U_xy) {
		return [U_xy[0] * U_xy[2],
				U_xy[0] * U_xy[1] * U_xy[2],
				U_xy[0] * U_xy[2] * U_xy[2] + this.g * U_xy[0] * U_xy[0] / 2];
	};
	// step function helpers
	this.U_half_h = function(i_lo, i_hi, j, dt) {
		var uhi = this.U(i_hi, j);
		var ulo = this.U(i_lo, j);
		var fhi = this.F(uhi);
		var flo = this.F(ulo);
		return [(uhi[0] + ulo[0]) * 0.5
				- (dt / (2 * this.dd)) * (fhi[0] - flo[0]),
				(uhi[1] + ulo[1]) * 0.5
				- (dt / (2 * this.dd)) * (fhi[1] - flo[1]),
				(uhi[2] + ulo[2]) * 0.5
				- (dt / (2 * this.dd)) * (fhi[2] - flo[2])];
	};
	this.U_half_v = function(i, j_lopez, j_hi, dt) {
		var uhi = this.U(i, j_hi);
		var ulo = this.U(i, j_lopez);
		var ghi = this.G(uhi);
		var glo = this.G(ulo);
		return [(uhi[0] + ulo[0]) * 0.5
				- (dt / (2 * this.dd)) * (ghi[0] - glo[0]),
				(uhi[1] + ulo[1]) * 0.5
				- (dt / (2 * this.dd)) * (ghi[1] - glo[1]),
				(uhi[2] + ulo[2]) * 0.5
				- (dt / (2 * this.dd)) * (ghi[2] - glo[2])];
	};
	// numerically step the system by time-diff dt
	// use Lax-Wendroff scheme above
	this.step = function(dt) {
		this.new_n = this.n;
		this.new_u = this.u;
		this.new_v = this.v;
		var uhhlo, uhhhi, uhvlo, uhvhi, new_U_xy;
		var fh, fl, gh, gl;
		var dtdd = dt / this.dd;
		for (var i = 0; i < this.QUANT; i++) {
			for (var j = 0; j < this.QUANT; j++) { // now in ass-polynomial time
				uhhlo = this.U_half_h(i - 1, i, j, dt);
				uhhhi = this.U_half_h(i, i + 1, j, dt);
				uhvlo = this.U_half_v(i, j - 1, j, dt);
				uhvhi = this.U_half_v(i, j, j + 1, dt);

				fh = this.F(uhhhi);
				fl = this.F(uhhlo);
				gh = this.G(uhvhi);
				gl = this.G(uhvlo);

				new_U_xy = [this.U(i,j)[0] - dtdd * (fh[0] - fl[0] + gh[0] - gl[0]),
							this.U(i,j)[1] - dtdd * (fh[1] - fl[1] + gh[1] - gl[1]),
							this.U(i,j)[2] - dtdd * (fh[2] - fl[2] + gh[2] - gl[2])];
				this.new_n[i][j] = new_U_xy[0];
				this.new_u[i][j] = new_U_xy[1] / new_U_xy[0];
				this.new_v[i][j] = new_U_xy[2] / new_U_xy[0];

				if (!isFinite((this.new_n[i][j]))) {
					console.log("bad new guy!!");
					console.trace();
					return;
				}
			}
		}
		this.n = this.new_n;
		this.u = this.new_u;
		this.v = this.new_v;
	};

	// it would be boring if we couldn't interact with this thing.
	// drip a drop height z at 0 < i,j < QUANT
	// guess it could be negative if you're really feelin it
	this.plip = function(i,j) {
		// console.log("plip",i,j);
		var drop = function(x, y) {
			// console.log("drop coords",(x-i),(y-j));
			return 30*((x - i) * (x - i) + (y - j) * (y - j));
		};
		var d;
		for (var p = 0; p < this.QUANT; p++) {
			for (var q = 0; q < this.QUANT; q++) {
				// console.log("dropping",drop(p,q),"at",p,q);
				d = drop(p,q);
				// console.log("old guy is",this.n[p][q],"adding",d);
				this.n[p][q] = this.n[p][q] + d;
				// console.log("new guy is",this.n[p][q]);
			}
		}
		return true; // for synchrony
	};

};

/* ********************************************************************************* */
/* visualization using canvas & simple colors                                        */
function CanvasBathtub($ctnr) {
	var Q = 500;

	this.width = 500;
	this.height = 500;
	this.canvas = document.createElement("canvas");
	this.canvas.setAttribute("id", "conway");
	this.canvas.setAttribute("width", this.width);
	this.canvas.setAttribute("height", this.height);
	$ctnr.append($(this.canvas));
	this.context = this.canvas.getContext("2d");
	this.context.imageSmoothingEnabled= false;
	this.id = undefined;
	var system = new SWE(Q, this.width);
	console.log(system);
	var self = this;

	var then = Date.now();
	var now;
	var dt;

	$(this.canvas).click(function(e) {
		var x = Math.floor((e.pageX-$(self.canvas).offset().left) / system.dd);
    	var y = Math.floor((e.pageY-$(self.canvas).offset().top) / system.dd);
    	// pause and drip a plip at x,y
    	self.cleanup();
    	if (system.plip(x,y))
    		self.loop();
	});

	this.animate = function() {
		console.log("ANIMATE DA FRAME",self.id);
		console.log(system);

		// render to canvas
		var b, c;
		for (var i = 0; i < Q; i++) {
			for (var j = 0; j < Q; j++) {
				if (isNaN((system.n[i][j]))) {
					console.log("bad!");
					console.trace();
					return;
				}
				b = (255 - Math.floor(system.n[i][j]) % 255).toString(16);
				c = "#" + ((b.length < 2) ? ("0" + b) : (b)) + ((b.length < 2) ? ("0" + b) : (b)) + "ff";
				// console.log("n-val",system.n[i][j],"produces color",c);
				self.context.fillStyle = c;
				self.context.fillRect(system.dd*i,system.dd*j,system.dd,system.dd);
			}
		}

		// timing for numerical step
		now = Date.now();
		system.step(now - then);
		then = now;

		self.id = requestAnimationFrame(self.animate);
	};
	this.loop = function() {
		this.id = requestAnimationFrame(self.animate, self.canvas);
	};
	this.cleanup = function() {
		cancelAnimationFrame(self.id);
	};
};























