// from http://en.wikipedia.org/wiki/Shallow_water_equations
// conservative form of shallow water equations:
// for 		n		- 	fluid column height at x,y,t
//			(u,v)	-	fluid horizontal velocity at x,y,t
//			g		-	acceleration due to gravity
//			H 		- 	mean water height
// if we write:
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
// then the system is:
//
//		dU     dF(U)     dG(U)     _
//		--  +  -----  +  -----  =  0
//		dt      dx        dy
//
// use Lax-Wendroff scheme per https://www.mathworks.com/moler/exm/chapters/water.pdf
// at any given moment, we have 2d array U_n, which has values at integer points (x,y)
// we want U_n+1. won't write out the method, but it's implemented below.

/* ********************************************************************************* */
// class SWE: interactive numerical shallow water system
// takes parameters above.
// reflective BCs: u = 0 at y = 0 or QUANT-1
//				   v = 0 at x = 0 or QUANT-1
// so our parameters are then:
//		QUANT		-	we're doing numerics. so we need a number of discrete points
//				    	to sample. our following parameters (except H) will thus be
//						QUANTxQUANT arrays.
//		LENGTH		-   the physical side length of the bathtub
// 		n_o			-	QUANTxQUANT array. initial value for column height.
//		u_o, v_o	-	same. IV for velocity. must conform to reflective BCs above.
function SWE(QUANT, LENGTH, n_o, u_o, v_o, g) {
	this.QUANT = QUANT || 100;
	this.LENGTH = LENGTH || 10;
	this.dd = this.LENGTH / this.QUANT;
	this.delegate; // will call updateHeightAt(x,y,height) on delegate whenever a height changes;
	// if n_o is not supplied, let's make a barely filled bathtub.
	if (n_o == undefined) {
		console.log("create n");
		n_o = [];
		for (var _i = 0; _i < QUANT; _i++) {
			n_o.push([]);
			for (var _j = 0; _j < QUANT; _j++) {
				n_o[_i].push(20.0);
			}
		}
	}
	if (u_o == undefined)
		u_o = n_o.map(function(a) { return a.map(function() { return 0.0; }) }); // backup is a correct-sized
	if (v_o == undefined)
		v_o = n_o.map(function(a) { return a.map(function() { return 0.0; }) }); // array of zeros
	// hold current state
	this.U = [];
	var h;
	for (var _ii = 0; _ii < QUANT; _ii++) {
		this.U.push([]);
		for (var _jj = 0; _jj < QUANT; _jj++) {
			h = n_o[_ii][_jj];
			this.U[_ii].push([h,
				h * u_o[_ii][_jj],
				h * v_o[_ii][_jj]]);
		}
	}
	this.g = g || 10; // meters per second^2
	// send current state
	this.heightAt = function(x,y) {
		return this.U[x][y][0];
	};
	// System
	this.F = function(U_xy) {
		var _hf = U_xy[0];
		var _uf = U_xy[1] / _hf;
		return [(U_xy[1]), 
		(U_xy[1] * _uf + this.g * _hf * _hf * 0.5),
		( U_xy[2] * _uf)];
	};
	this.G = function(U_xy) {
		var _hg = U_xy[0];
		var _vg = U_xy[2] / _hg;
		return [(U_xy[2]),
		(U_xy[1] * _vg),
		(U_xy[2] * _vg + this.g * _hg * _hg * 0.5)];
	};
	// step function helpers
	function addthree(a,b) {
		return [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
	}
	function subthree(a,b) {
		return [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
	}
	function scalethree(a,vec) {
		return [a*vec[0],a*vec[1],a*vec[2]];
	}
	// numerically step the system by time-diff dt (seconds)
	// use Lax-Wendroff scheme above
	// this is a meaty method yall! now coming at you in ass-polynomial time.
	// whence the NaNs??
	this.step = function(dt) {
		// have U^{n}_{i,j} of size QxQ stored in this.U
		// that gives values here:
		// |===|
		// | x |
		// |===|
		// remember BCs: u = 0 at y = 0 or QUANT-1
		//				 v = 0 at x = 0 or QUANT-1


		/*******************/
		/* STEP 1 Lax shiz */
		/*******************/
		var step_factor = dt * 0.5 / this.dd;
		// this will be U^{n+.5}_{i+.5,j} of size Q+1xQ - the horizontal midpoints
		// can see the size change cuz it gives values here:
		// |===|
		// x   x
		// |===|
		// careful with cases outside of U's dims:
		// ih = 0 (i = -.5 in U coords) and
		// ih = Q (i = Q - .5 in U coords)
		// handle these by substituting U values for i = 0 and i = Q - 1
		// also make sure 3rd component is 0 at these points
		// for vertical edges, make sure U^{n+.5}_{i+.5,0} and U^{n+.5}_{i+.5,Q-1}
		// have second component equal to zero.
		var U_lax_horizontal = [];
		var uhlo, uhhi; // store points to be interpolated
		for (var ih = 0; ih < this.QUANT + 1; ih++) {
			// ih is i + 1/2 in original coords
			U_lax_horizontal.push([]);
			for (var jh = 0; jh < this.QUANT; jh++) {
				// jh is j in original coords
				// interpolate U[ih - 1][jh] and U[ih][jh]
				
				// edge cases:
				if (ih == 0) {
					uhlo = this.U[ih][jh];
					uhhi = uhlo;
				} else if (ih == this.QUANT) {
					uhlo = this.U[ih - 1][jh];
					uhhi = uhlo;
				} else {
					uhlo = this.U[ih - 1][jh];
					uhhi = this.U[ih][jh];
				}
				// Lax interpolate
				U_lax_horizontal[ih][jh] = subthree(scalethree(0.5, addthree(uhhi, uhlo)),
					scalethree(step_factor,
						subthree(this.F(uhhi), this.F(uhlo))));
			}
		}

		// this will be U^{n+.5}_{i,j+.5} of size QxQ+1 - the vertical midpoints
		// size change cuz it gives values here:
		// |=x=|
		// |   |
		// |=x=|
		// careful with edge cases:
		// jv = 0 (j = -1/2 in this.U coords) and
		// jv = Q+1 (j = Q + .5 in U coords)
		var U_lax_vertical = [];
		var uvlo, uvhi; // store points to be interpolated
		for (var iv = 0; iv < this.QUANT; iv++) {
			// iv is i in original coords
			U_lax_vertical.push([]);
			for (var jv = 0; jv < this.QUANT + 1; jv++) {
				// jv is j + 1/2 in original coords
				
				// edges:
				if (jv == 0) {
					uvhi = this.U[iv][jv];
					uvlo = uvhi;
				} else if (jv == this.QUANT) {
					uvlo = this.U[iv][jv - 1];
					uvhi = uvlo;
				} else {
					uvlo = this.U[iv][jv - 1];
					uvhi = this.U[iv][jv];
				}
				// Lax interpolate
				U_lax_vertical[iv][jv] = subthree((scalethree(0.5, addthree(uvhi, uvlo))),
					scalethree(step_factor,
						subthree(this.G(uvhi), this.G(uvlo))));
			}
		}

		/*****************/
		/* Step 2 U_next */
		/*****************/
		// we're now generating U^{n+1}_{i,j} QxQ, the new value of U
		var next_factor = step_factor * 2;
		var U_next = [];
		for (var i = 0; i < this.QUANT; i++) {
			U_next.push([]);
			for (var j = 0; j < this.QUANT; j++) {
				// going to work with:
				// U^{n+.5}_{i+.5,j} and U^{n+.5}_{i-.5,j}
				//		find with [i+1][j] and [i][j], respectively
				// U^{n+.5}_{i,j+.5} and U^{n+.5}_{i,j-.5}
				//		find with [i][j+j1] and [i][j], respectively

				// functional programming has ruined me, dear readers
				U_next[i].push(subthree(this.U[i][j],
					scalethree(next_factor,
						addthree(subthree(this.F(U_lax_horizontal[i+1][j]),
							this.F(U_lax_horizontal[i][j])),
						subthree(this.G(U_lax_vertical[i][j+1]),
							this.G(U_lax_vertical[i][j]))))));
				if (this.delegate && U_next[i][j][0] != this.U[i][j][0])
					this.delegate.updateHeightAt(i,j,U_next[i][j][0]);
			}
		}
		this.U = U_next; // ay
		return true;
	};

	// plip helpers
	this.addtoU = function(a,i,j) {
		// got to break it down
		this.U[i][j][1] = this.U[i][j][1] / this.U[i][j][0];
		this.U[i][j][2] = this.U[i][j][2] / this.U[i][j][0];
		this.U[i][j][0] += a;
		this.U[i][j][1] *= this.U[i][j][0];
		this.U[i][j][2] *= this.U[i][j][0];
		for (var p = 0; p < 3; p++) {
			console.assert(isFinite(this.U[i][j][p]));
		}
	};
	var drop = function(x, y, i, j) {
		return 100 * Math.exp(-0.008 * ((x - i) * (x - i) + (y - j) * (y - j)));
	};
	// it would be boring if we couldn't interact with this thing.
	// drip a drop at 0 < i,j < QUANT
	this.plip = function(i,j) {
		console.log("PLIP",i,j);
		var d;
		for (var p = 0; p < this.QUANT; p++) {
			for (var q = 0; q < this.QUANT; q++) {
				d = drop(p,q,i,j);
				this.addtoU(d,p,q);
			}
		}
		return true;
	};

	this.randomize = function() {
		this.U = this.U.map(function (a) { return a.map(function(b) { return b.map(function() { return Math.random() * 40; })})});
	}

};

/* ********************************************************************************* */
/* visualization using canvas & simple colors                                        */
function CanvasBathtub($ctnr) {
	var Q = 150;

	this.width = 750;
	this.height = 750;
	this.canvas = document.createElement("canvas");
	this.canvas.setAttribute("id", "conway");
	this.canvas.setAttribute("width", this.width);
	this.canvas.setAttribute("height", this.height);
	$ctnr.append($(this.canvas));
	this.context = this.canvas.getContext("2d");
	this.context.imageSmoothingEnabled = false;
	this.running = true;
	var system = new SWE(Q, this.width);
	console.log(system);
	var self = this;

	var then = Date.now();
	var now;
	var dt;

	$(this.canvas).click(function(e) {
		var x = Math.floor(Math.floor((e.pageX-$(self.canvas).offset().left)) / (system.dd));
		var y = Math.floor(Math.floor((e.pageY-$(self.canvas).offset().top)) / (system.dd));
    	// pause and drip a plip at x,y
    	self.running = false;
    	system.plip(x,y);
    	then = Date.now() - 1;
    	self.update();
    	self.running = true;
    	self.loop();
    });

	// update to current tub.
	this.update = function() {
		console.log("refresh!");
		var height;
		for (var x = 0; x < Q; x++) {
			for (var y = 0; y < Q; y++) {
				height = system.heightAt(x,y);
				var col = (255 - (Math.round(height) % 256));
				this.context.fillStyle = "rgb("+ col + ","+ col +",255)";
				this.context.fillRect(system.dd * x, system.dd * y, system.dd, system.dd);
			}
		}
	};

	// be a bathtub delegate
	system.delegate = this;
	this.updateHeightAt = function(x,y,height) {
		var o = 
		console.assert(isFinite(height));
		var col = (255 - (Math.round(height) % 256));
		this.context.fillStyle = "rgb("+ col + ","+ col +",255)";
		if (!isFinite(height)) {
			console.log("bimmer at",x,y,height);
			console.log("color was",this.context.fillStyle);
			this.context.fillStyle = "#000";
		}
		this.context.fillRect(system.dd * x, system.dd * y, system.dd, system.dd);
	};

	this.loop = function() {
		console.log("loop")
		now = Date.now();
		dt = now - then;
		system.step(dt / 10000);
		then = now;
		if (self.running) {
			window.requestAnimationFrame(self.loop);
		}
	};

	this.cleanup = function() { 
		this.running = false;
	};
	this.update();
};
























function Conway($container) {
	var self = this; // love JS!
	this.state = [];
	this.width = Math.round($container.width());
	this.height = Math.round($container.height());
	this.loopid = null;

	// init empty board
	for (var x = 0; x < this.width; x++) {
		this.state.push([]);
		for (var y = 0; y < this.height; y++) {
			this.state[x].push(false);
		}
	};

	this.$container = $container;
	// init canvas
	this.canvas = document.createElement("canvas");
	this.canvas.setAttribute("id", "conway");
	this.canvas.setAttribute("width", this.width);
	this.canvas.setAttribute("height", this.height);
	this.$container.append($(this.canvas));
	this.context = this.canvas.getContext("2d");
	this.context.imageSmoothingEnabled= false;

	// get image data from current state
	this.toImage = function() {
		var imgd = this.context.createImageData(this.width, this.height);
		for (var x = 0; x < this.width; x++) {
			for (var y = 0; y < this.height; y++) {
				if (this.state[x][y]) {
					imgd.data[4 * x + 4 * this.width * y] = 0;
					imgd.data[4 * x + 4 * this.width * y + 1] = 0;
					imgd.data[4 * x + 4 * this.width * y + 2] = 0;
					imgd.data[4 * x + 4 * this.width * y + 3] = 255;
				} else {
					imgd.data[4 * x + 4 * this.width * y + 1] = 255;
					imgd.data[4 * x + 4 * this.width * y + 2] = 255;
					imgd.data[4 * x + 4 * this.width * y + 3] = 255;
					imgd.data[4 * x + 4 * this.width * y] = 255;
				}
			}
		}
		return imgd;
	};

	// insert a width x height canvas into container
	// seed the game
	this.loop = function() {
		var self = this;
		this.loopid = setInterval(function() {
			console.log("looping");
			self.context.putImageData(self.toImage(), 0, 0);
			self.next();
		}, 100);
	};

	// cleanup tasks
	this.destroy = function() {
		clearInterval(self.loopid);
	};


	// fill the state with falses
	this.clear = function() {
		for (var x = 0; x < this.width; x++) {
			for (var y = 0; y < this.height; y++) {
				this.state[x][y] = false;
			}
		}
	}; // end clear

	// refill the state randomly
	this.randomSeed = function() {
		for (var x = 0; x < this.width; x++) {
			for (var y = 0; y < this.height; y++) {
						this.state[x][y] = Math.random() < .5;; // increment if the cell is true
			}
		}
	}; // end randomSeed
	this.randomSeed();

	// count all of the living cells
	this.count = function() {
		var count = 0;
		for (var x = 0; x < this.width; x++) {
			for (var y = 0; y < this.height; y++) {
				if (this.state[x][y]) {
						count++; // increment if the cell is true
					}
				}
			}
			return count;
	}; // end count

	// transition to the next state
	// if a cell is alive, it stays alive if it has 2 or 3 neighbors
	// a dead cell comes to life with exactly 3 neighbors
	this.next = function() {
		var next_state = new Array(this.width);
		for (var x = 0; x < this.width; x++) {
			next_state[x] = new Array(this.height);
		}

		for (var x = 0; x < this.width; x++) {
			for (var y = 0; y < this.height; y++) {
				next_state[x][y] = this.lives(x, y);
			}
		}
		this.state = next_state;
		return this.state;
	}; // end next

	// see if the cell at x, y will live to the next generation
	this.lives = function(x, y) {
		var neighbors = this.neighbors(x, y);
		if (this.state[x][y]) {
			return (neighbors > 1 && neighbors < 4);
		} else {
			return (neighbors == 3);
		}
	}; // end lives

	// count a cell's living neighbors
	this.neighbors = function(x, y) {
		var n = 0;
		if (x > 0 && x < this.width - 1 && y > 0 && y < this.height) {
			for (var i = -1; i <= 1; i++) {
				for (var p = -1; p <= 1; p++) {
					if (!(i == 0 && p == 0) && this.state[x + i][y + p]) n++;
				}
			}
		} else if (x == 0 && y == 0) {
			if (this.state[x][y + 1]) n++;
			if (this.state[x + 1][y]) n++;
			if (this.state[x + 1][y + 1]) n++;
		} else if (x == 0 && y == this.height - 1) {
			if (this.state[x][y - 1]) n++;
			if (this.state[x + 1][y - 1]) n++;
			if (this.state[x + 1][y]) n++;
		} else if (x == this.width - 1 && y == 0) {
			if (this.state[x - 1][y]) n++;
			if (this.state[x - 1][y + 1]) n++;
			if (this.state[x][y + 1]) n++;
		} else if (x == this.width - 1 && y == this.height - 1) {
			if (this.state[x][y - 1]) n++;
			if (this.state[x - 1][y]) n++;
			if (this.state[x - 1][y - 1]) n++;
		} else if (x == 0) {
			if (this.state[x][y + 1]) n++;
			if (this.state[x][y - 1]) n++;
			if (this.state[x + 1][y + 1]) n++;
			if (this.state[x + 1][y - 1]) n++;
			if (this.state[x + 1][y]) n++;
		} else if (x == this.width - 1) {
			if (this.state[x][y + 1]) n++;
			if (this.state[x][y - 1]) n++;
			if (this.state[x - 1][y + 1]) n++;
			if (this.state[x - 1][y - 1]) n++;
			if (this.state[x - 1][y]) n++;
		} else if (y == 0) {
			if (this.state[x + 1][y]) n++;
			if (this.state[x - 1][y]) n++;
			if (this.state[x + 1][y + 1]) n++;
			if (this.state[x - 1][y + 1]) n++;
			if (this.state[x][y + 1]) n++;
		} else if (y == this.height - 1) {
			if (this.state[x + 1][y]) n++;
			if (this.state[x - 1][y]) n++;
			if (this.state[x + 1][y - 1]) n++;
			if (this.state[x - 1][y - 1]) n++;
			if (this.state[x][y - 1]) n++;
		}
		return n;
	}; // end neighbors

	this.toString = function() {
		var out = "";
		for (var y = 0; y < this.height; y++) {
			for (var x = 0; x < this.width; x++) {
				if (this.state[x][y]) {
					out = out + "+ ";
				} else {
					out = out + "  ";
				}
			}
			out = out + "\n";
		}
		return out;
	}; // /tostring
	// drop a glider randomly on the grid
	this.glider = function() {
			console.log("glider: " + this);
			var xcenter = Math.floor(Math.random() * (this.width * 0.8 + 1) + this.width * 0.1);
			var ycenter = Math.floor(Math.random() * (this.height * 0.8 + 1) + this.height * 0.1);
			this.state[xcenter][ycenter + 1] = true;
			this.state[xcenter - 1][ycenter] = true;
			this.state[xcenter - 1][ycenter - 1] = true;
			this.state[xcenter][ycenter - 1] = true;
			this.state[xcenter + 1][ycenter - 1] = true;

	}; // /glider
	// drop a glider gun randomly on the grid
	this.gliderGun = function() {
			var xcenter = Math.floor(Math.random() * (this.width * 0.8 + 1) + this.width * 0.1);
			var ycenter = Math.floor(Math.random() * (this.height * 0.8 + 1) + this.height * 0.1);
			this.state[xcenter + 6][ycenter + -4] = true;
			this.state[xcenter + 4][ycenter + -3] = true;
			this.state[xcenter + 6][ycenter + -3] = true;
			this.state[xcenter + -6][ycenter + -2] = true;
			this.state[xcenter + -5][ycenter + -2] = true;
			this.state[xcenter + 2][ycenter + -2] = true;
			this.state[xcenter + 3][ycenter + -2] = true;
			this.state[xcenter + 16][ycenter + -2] = true;
			this.state[xcenter + 17][ycenter + -2] = true;
			this.state[xcenter + -7][ycenter + -1] = true;
			this.state[xcenter + -3][ycenter + -1] = true;
			this.state[xcenter + 2][ycenter + -1] = true;
			this.state[xcenter + 3][ycenter + -1] = true;
			this.state[xcenter + 16][ycenter + -1] = true;
			this.state[xcenter + 17][ycenter + -1] = true;
			this.state[xcenter + -18][ycenter + 0] = true;
			this.state[xcenter + -17][ycenter + 0] = true;
			this.state[xcenter + -8][ycenter + 0] = true;
			this.state[xcenter + -2][ycenter + 0] = true;
			this.state[xcenter + 2][ycenter + 0] = true;
			this.state[xcenter + 3][ycenter + 0] = true;
			this.state[xcenter + -18][ycenter + 1] = true;
			this.state[xcenter + -17][ycenter + 1] = true;
			this.state[xcenter + -8][ycenter + 1] = true;
			this.state[xcenter + -4][ycenter + 1] = true;
			this.state[xcenter + -2][ycenter + 1] = true;
			this.state[xcenter + -1][ycenter + 1] = true;
			this.state[xcenter + 4][ycenter + 1] = true;
			this.state[xcenter + 6][ycenter + 1] = true;
			this.state[xcenter + -8][ycenter + 2] = true;
			this.state[xcenter + -2][ycenter + 2] = true;
			this.state[xcenter + 6][ycenter + 2] = true;
			this.state[xcenter + -7][ycenter + 3] = true;
			this.state[xcenter + -3][ycenter + 3] = true;
			this.state[xcenter + -6][ycenter + 4] = true;
			this.state[xcenter + -5][ycenter + 4] = true;
	};
	// drop puffer randomly on the grid
	this.puffer = function() {
			var xcenter = Math.floor(Math.random() * (this.width * 0.7 + 1) + this.width * 0.15);
			var ycenter = Math.floor(Math.random() * (this.height * 0.8 + 1) + this.height * 0.1);
			this.state[xcenter + 0][ycenter + 43] = true;
			this.state[xcenter + 1][ycenter + 41] = true;
			this.state[xcenter + 1][ycenter + 43] = true;
			this.state[xcenter + 2][ycenter + 42] = true;
			this.state[xcenter + 2][ycenter + 43] = true;
			this.state[xcenter + 3][ycenter + 84] = true;
			this.state[xcenter + 4][ycenter + 82] = true;
			this.state[xcenter + 4][ycenter + 83] = true;
			this.state[xcenter + 5][ycenter + 83] = true;
			this.state[xcenter + 5][ycenter + 84] = true;
			this.state[xcenter + 10][ycenter + 58] = true;
			this.state[xcenter + 10][ycenter + 59] = true;
			this.state[xcenter + 11][ycenter + 57] = true;
			this.state[xcenter + 11][ycenter + 59] = true;
			this.state[xcenter + 12][ycenter + 46] = true;
			this.state[xcenter + 12][ycenter + 47] = true;
			this.state[xcenter + 12][ycenter + 57] = true;
			this.state[xcenter + 12][ycenter + 58] = true;
			this.state[xcenter + 12][ycenter + 68] = true;
			this.state[xcenter + 12][ycenter + 69] = true;
			this.state[xcenter + 13][ycenter + 46] = true;
			this.state[xcenter + 13][ycenter + 47] = true;
			this.state[xcenter + 13][ycenter + 68] = true;
			this.state[xcenter + 13][ycenter + 70] = true;
			this.state[xcenter + 14][ycenter + 69] = true;
			this.state[xcenter + 14][ycenter + 70] = true;
			this.state[xcenter + 14][ycenter + 80] = true;
			this.state[xcenter + 14][ycenter + 81] = true;
			this.state[xcenter + 15][ycenter + 28] = true;
			this.state[xcenter + 15][ycenter + 80] = true;
			this.state[xcenter + 15][ycenter + 81] = true;
			this.state[xcenter + 16][ycenter + 26] = true;
			this.state[xcenter + 16][ycenter + 28] = true;
			this.state[xcenter + 17][ycenter + 15] = true;
			this.state[xcenter + 17][ycenter + 17] = true;
			this.state[xcenter + 17][ycenter + 27] = true;
			this.state[xcenter + 17][ycenter + 28] = true;
			this.state[xcenter + 18][ycenter + 14] = true;
			this.state[xcenter + 18][ycenter + 99] = true;
			this.state[xcenter + 19][ycenter + 14] = true;
			this.state[xcenter + 19][ycenter + 18] = true;
			this.state[xcenter + 19][ycenter + 97] = true;
			this.state[xcenter + 19][ycenter + 98] = true;
			this.state[xcenter + 20][ycenter + 14] = true;
			this.state[xcenter + 20][ycenter + 98] = true;
			this.state[xcenter + 20][ycenter + 99] = true;
			this.state[xcenter + 21][ycenter + 14] = true;
			this.state[xcenter + 21][ycenter + 17] = true;
			this.state[xcenter + 21][ycenter + 110] = true;
			this.state[xcenter + 21][ycenter + 112] = true;
			this.state[xcenter + 22][ycenter + 14] = true;
			this.state[xcenter + 22][ycenter + 15] = true;
			this.state[xcenter + 22][ycenter + 16] = true;
			this.state[xcenter + 22][ycenter + 113] = true;
			this.state[xcenter + 23][ycenter + 109] = true;
			this.state[xcenter + 23][ycenter + 113] = true;
			this.state[xcenter + 24][ycenter + 113] = true;
			this.state[xcenter + 25][ycenter + 18] = true;
			this.state[xcenter + 25][ycenter + 110] = true;
			this.state[xcenter + 25][ycenter + 113] = true;
			this.state[xcenter + 26][ycenter + 16] = true;
			this.state[xcenter + 26][ycenter + 18] = true;
			this.state[xcenter + 26][ycenter + 22] = true;
			this.state[xcenter + 26][ycenter + 111] = true;
			this.state[xcenter + 26][ycenter + 112] = true;
			this.state[xcenter + 26][ycenter + 113] = true;
			this.state[xcenter + 27][ycenter + 17] = true;
			this.state[xcenter + 27][ycenter + 18] = true;
			this.state[xcenter + 27][ycenter + 22] = true;
			this.state[xcenter + 28][ycenter + 43] = true;
			this.state[xcenter + 28][ycenter + 109] = true;
			this.state[xcenter + 29][ycenter + 43] = true;
			this.state[xcenter + 29][ycenter + 44] = true;
			this.state[xcenter + 29][ycenter + 107] = true;
			this.state[xcenter + 29][ycenter + 108] = true;
			this.state[xcenter + 30][ycenter + 13] = true;
			this.state[xcenter + 30][ycenter + 43] = true;
			this.state[xcenter + 30][ycenter + 44] = true;
			this.state[xcenter + 30][ycenter + 45] = true;
			this.state[xcenter + 30][ycenter + 84] = true;
			this.state[xcenter + 30][ycenter + 108] = true;
			this.state[xcenter + 30][ycenter + 109] = true;
			this.state[xcenter + 31][ycenter + 11] = true;
			this.state[xcenter + 31][ycenter + 13] = true;
			this.state[xcenter + 31][ycenter + 19] = true;
			this.state[xcenter + 31][ycenter + 20] = true;
			this.state[xcenter + 31][ycenter + 24] = true;
			this.state[xcenter + 31][ycenter + 25] = true;
			this.state[xcenter + 31][ycenter + 37] = true;
			this.state[xcenter + 31][ycenter + 45] = true;
			this.state[xcenter + 31][ycenter + 46] = true;
			this.state[xcenter + 31][ycenter + 83] = true;
			this.state[xcenter + 31][ycenter + 84] = true;
			this.state[xcenter + 32][ycenter + 12] = true;
			this.state[xcenter + 32][ycenter + 13] = true;
			this.state[xcenter + 32][ycenter + 36] = true;
			this.state[xcenter + 32][ycenter + 37] = true;
			this.state[xcenter + 32][ycenter + 38] = true;
			this.state[xcenter + 32][ycenter + 40] = true;
			this.state[xcenter + 32][ycenter + 44] = true;
			this.state[xcenter + 32][ycenter + 45] = true;
			this.state[xcenter + 32][ycenter + 82] = true;
			this.state[xcenter + 32][ycenter + 83] = true;
			this.state[xcenter + 32][ycenter + 84] = true;
			this.state[xcenter + 33][ycenter + 35] = true;
			this.state[xcenter + 33][ycenter + 36] = true;
			this.state[xcenter + 33][ycenter + 38] = true;
			this.state[xcenter + 33][ycenter + 39] = true;
			this.state[xcenter + 33][ycenter + 40] = true;
			this.state[xcenter + 33][ycenter + 41] = true;
			this.state[xcenter + 33][ycenter + 43] = true;
			this.state[xcenter + 33][ycenter + 44] = true;
			this.state[xcenter + 33][ycenter + 57] = true;
			this.state[xcenter + 33][ycenter + 58] = true;
			this.state[xcenter + 33][ycenter + 81] = true;
			this.state[xcenter + 33][ycenter + 82] = true;
			this.state[xcenter + 33][ycenter + 90] = true;
			this.state[xcenter + 33][ycenter + 114] = true;
			this.state[xcenter + 34][ycenter + 22] = true;
			this.state[xcenter + 34][ycenter + 29] = true;
			this.state[xcenter + 34][ycenter + 30] = true;
			this.state[xcenter + 34][ycenter + 36] = true;
			this.state[xcenter + 34][ycenter + 38] = true;
			this.state[xcenter + 34][ycenter + 39] = true;
			this.state[xcenter + 34][ycenter + 41] = true;
			this.state[xcenter + 34][ycenter + 42] = true;
			this.state[xcenter + 34][ycenter + 56] = true;
			this.state[xcenter + 34][ycenter + 59] = true;
			this.state[xcenter + 34][ycenter + 82] = true;
			this.state[xcenter + 34][ycenter + 83] = true;
			this.state[xcenter + 34][ycenter + 87] = true;
			this.state[xcenter + 34][ycenter + 89] = true;
			this.state[xcenter + 34][ycenter + 90] = true;
			this.state[xcenter + 34][ycenter + 91] = true;
			this.state[xcenter + 34][ycenter + 112] = true;
			this.state[xcenter + 34][ycenter + 113] = true;
			this.state[xcenter + 35][ycenter + 8] = true;
			this.state[xcenter + 35][ycenter + 21] = true;
			this.state[xcenter + 35][ycenter + 23] = true;
			this.state[xcenter + 35][ycenter + 28] = true;
			this.state[xcenter + 35][ycenter + 29] = true;
			this.state[xcenter + 35][ycenter + 30] = true;
			this.state[xcenter + 35][ycenter + 37] = true;
			this.state[xcenter + 35][ycenter + 38] = true;
			this.state[xcenter + 35][ycenter + 39] = true;
			this.state[xcenter + 35][ycenter + 41] = true;
			this.state[xcenter + 35][ycenter + 46] = true;
			this.state[xcenter + 35][ycenter + 57] = true;
			this.state[xcenter + 35][ycenter + 58] = true;
			this.state[xcenter + 35][ycenter + 69] = true;
			this.state[xcenter + 35][ycenter + 70] = true;
			this.state[xcenter + 35][ycenter + 83] = true;
			this.state[xcenter + 35][ycenter + 84] = true;
			this.state[xcenter + 35][ycenter + 86] = true;
			this.state[xcenter + 35][ycenter + 87] = true;
			this.state[xcenter + 35][ycenter + 88] = true;
			this.state[xcenter + 35][ycenter + 89] = true;
			this.state[xcenter + 35][ycenter + 91] = true;
			this.state[xcenter + 35][ycenter + 92] = true;
			this.state[xcenter + 35][ycenter + 105] = true;
			this.state[xcenter + 35][ycenter + 113] = true;
			this.state[xcenter + 35][ycenter + 114] = true;
			this.state[xcenter + 36][ycenter + 6] = true;
			this.state[xcenter + 36][ycenter + 8] = true;
			this.state[xcenter + 36][ycenter + 21] = true;
			this.state[xcenter + 36][ycenter + 23] = true;
			this.state[xcenter + 36][ycenter + 28] = true;
			this.state[xcenter + 36][ycenter + 29] = true;
			this.state[xcenter + 36][ycenter + 31] = true;
			this.state[xcenter + 36][ycenter + 38] = true;
			this.state[xcenter + 36][ycenter + 39] = true;
			this.state[xcenter + 36][ycenter + 45] = true;
			this.state[xcenter + 36][ycenter + 47] = true;
			this.state[xcenter + 36][ycenter + 68] = true;
			this.state[xcenter + 36][ycenter + 71] = true;
			this.state[xcenter + 36][ycenter + 85] = true;
			this.state[xcenter + 36][ycenter + 86] = true;
			this.state[xcenter + 36][ycenter + 88] = true;
			this.state[xcenter + 36][ycenter + 89] = true;
			this.state[xcenter + 36][ycenter + 91] = true;
			this.state[xcenter + 36][ycenter + 97] = true;
			this.state[xcenter + 36][ycenter + 98] = true;
			this.state[xcenter + 36][ycenter + 104] = true;
			this.state[xcenter + 36][ycenter + 105] = true;
			this.state[xcenter + 36][ycenter + 106] = true;
			this.state[xcenter + 37][ycenter + 7] = true;
			this.state[xcenter + 37][ycenter + 8] = true;
			this.state[xcenter + 37][ycenter + 22] = true;
			this.state[xcenter + 37][ycenter + 29] = true;
			this.state[xcenter + 37][ycenter + 30] = true;
			this.state[xcenter + 37][ycenter + 31] = true;
			this.state[xcenter + 37][ycenter + 45] = true;
			this.state[xcenter + 37][ycenter + 47] = true;
			this.state[xcenter + 37][ycenter + 48] = true;
			this.state[xcenter + 37][ycenter + 69] = true;
			this.state[xcenter + 37][ycenter + 70] = true;
			this.state[xcenter + 37][ycenter + 81] = true;
			this.state[xcenter + 37][ycenter + 86] = true;
			this.state[xcenter + 37][ycenter + 88] = true;
			this.state[xcenter + 37][ycenter + 89] = true;
			this.state[xcenter + 37][ycenter + 90] = true;
			this.state[xcenter + 37][ycenter + 97] = true;
			this.state[xcenter + 37][ycenter + 98] = true;
			this.state[xcenter + 37][ycenter + 99] = true;
			this.state[xcenter + 37][ycenter + 104] = true;
			this.state[xcenter + 37][ycenter + 105] = true;
			this.state[xcenter + 37][ycenter + 106] = true;
			this.state[xcenter + 38][ycenter + 22] = true;
			this.state[xcenter + 38][ycenter + 30] = true;
			this.state[xcenter + 38][ycenter + 46] = true;
			this.state[xcenter + 38][ycenter + 47] = true;
			this.state[xcenter + 38][ycenter + 80] = true;
			this.state[xcenter + 38][ycenter + 82] = true;
			this.state[xcenter + 38][ycenter + 88] = true;
			this.state[xcenter + 38][ycenter + 89] = true;
			this.state[xcenter + 38][ycenter + 96] = true;
			this.state[xcenter + 38][ycenter + 98] = true;
			this.state[xcenter + 38][ycenter + 99] = true;
			this.state[xcenter + 38][ycenter + 102] = true;
			this.state[xcenter + 38][ycenter + 104] = true;
			this.state[xcenter + 38][ycenter + 106] = true;
			this.state[xcenter + 38][ycenter + 108] = true;
			this.state[xcenter + 38][ycenter + 119] = true;
			this.state[xcenter + 39][ycenter + 21] = true;
			this.state[xcenter + 39][ycenter + 23] = true;
			this.state[xcenter + 39][ycenter + 79] = true;
			this.state[xcenter + 39][ycenter + 80] = true;
			this.state[xcenter + 39][ycenter + 82] = true;
			this.state[xcenter + 39][ycenter + 96] = true;
			this.state[xcenter + 39][ycenter + 97] = true;
			this.state[xcenter + 39][ycenter + 98] = true;
			this.state[xcenter + 39][ycenter + 103] = true;
			this.state[xcenter + 39][ycenter + 107] = true;
			this.state[xcenter + 39][ycenter + 117] = true;
			this.state[xcenter + 39][ycenter + 118] = true;
			this.state[xcenter + 40][ycenter + 20] = true;
			this.state[xcenter + 40][ycenter + 24] = true;
			this.state[xcenter + 40][ycenter + 58] = true;
			this.state[xcenter + 40][ycenter + 59] = true;
			this.state[xcenter + 40][ycenter + 80] = true;
			this.state[xcenter + 40][ycenter + 81] = true;
			this.state[xcenter + 40][ycenter + 97] = true;
			this.state[xcenter + 40][ycenter + 103] = true;
			this.state[xcenter + 40][ycenter + 107] = true;
			this.state[xcenter + 40][ycenter + 118] = true;
			this.state[xcenter + 40][ycenter + 119] = true;
			this.state[xcenter + 41][ycenter + 21] = true;
			this.state[xcenter + 41][ycenter + 23] = true;
			this.state[xcenter + 41][ycenter + 41] = true;
			this.state[xcenter + 41][ycenter + 47] = true;
			this.state[xcenter + 41][ycenter + 57] = true;
			this.state[xcenter + 41][ycenter + 59] = true;
			this.state[xcenter + 42][ycenter + 1] = true;
			this.state[xcenter + 42][ycenter + 3] = true;
			this.state[xcenter + 42][ycenter + 22] = true;
			this.state[xcenter + 42][ycenter + 41] = true;
			this.state[xcenter + 42][ycenter + 42] = true;
			this.state[xcenter + 42][ycenter + 47] = true;
			this.state[xcenter + 42][ycenter + 57] = true;
			this.state[xcenter + 42][ycenter + 58] = true;
			this.state[xcenter + 42][ycenter + 68] = true;
			this.state[xcenter + 42][ycenter + 69] = true;
			this.state[xcenter + 43][ycenter + 0] = true;
			this.state[xcenter + 43][ycenter + 40] = true;
			this.state[xcenter + 43][ycenter + 42] = true;
			this.state[xcenter + 43][ycenter + 47] = true;
			this.state[xcenter + 43][ycenter + 68] = true;
			this.state[xcenter + 43][ycenter + 70] = true;
			this.state[xcenter + 43][ycenter + 80] = true;
			this.state[xcenter + 43][ycenter + 86] = true;
			this.state[xcenter + 44][ycenter + 0] = true;
			this.state[xcenter + 44][ycenter + 19] = true;
			this.state[xcenter + 44][ycenter + 20] = true;
			this.state[xcenter + 44][ycenter + 24] = true;
			this.state[xcenter + 44][ycenter + 25] = true;
			this.state[xcenter + 44][ycenter + 69] = true;
			this.state[xcenter + 44][ycenter + 70] = true;
			this.state[xcenter + 44][ycenter + 80] = true;
			this.state[xcenter + 44][ycenter + 85] = true;
			this.state[xcenter + 44][ycenter + 86] = true;
			this.state[xcenter + 44][ycenter + 104] = true;
			this.state[xcenter + 44][ycenter + 105] = true;
			this.state[xcenter + 44][ycenter + 106] = true;
			this.state[xcenter + 45][ycenter + 0] = true;
			this.state[xcenter + 45][ycenter + 3] = true;
			this.state[xcenter + 45][ycenter + 10] = true;
			this.state[xcenter + 45][ycenter + 11] = true;
			this.state[xcenter + 45][ycenter + 12] = true;
			this.state[xcenter + 45][ycenter + 19] = true;
			this.state[xcenter + 45][ycenter + 20] = true;
			this.state[xcenter + 45][ycenter + 21] = true;
			this.state[xcenter + 45][ycenter + 23] = true;
			this.state[xcenter + 45][ycenter + 24] = true;
			this.state[xcenter + 45][ycenter + 25] = true;
			this.state[xcenter + 45][ycenter + 46] = true;
			this.state[xcenter + 45][ycenter + 48] = true;
			this.state[xcenter + 45][ycenter + 80] = true;
			this.state[xcenter + 45][ycenter + 85] = true;
			this.state[xcenter + 45][ycenter + 87] = true;
			this.state[xcenter + 45][ycenter + 104] = true;
			this.state[xcenter + 45][ycenter + 105] = true;
			this.state[xcenter + 45][ycenter + 106] = true;
			this.state[xcenter + 46][ycenter + 0] = true;
			this.state[xcenter + 46][ycenter + 1] = true;
			this.state[xcenter + 46][ycenter + 2] = true;
			this.state[xcenter + 46][ycenter + 10] = true;
			this.state[xcenter + 46][ycenter + 11] = true;
			this.state[xcenter + 46][ycenter + 12] = true;
			this.state[xcenter + 46][ycenter + 18] = true;
			this.state[xcenter + 46][ycenter + 20] = true;
			this.state[xcenter + 46][ycenter + 21] = true;
			this.state[xcenter + 46][ycenter + 23] = true;
			this.state[xcenter + 46][ycenter + 24] = true;
			this.state[xcenter + 46][ycenter + 26] = true;
			this.state[xcenter + 46][ycenter + 45] = true;
			this.state[xcenter + 46][ycenter + 46] = true;
			this.state[xcenter + 46][ycenter + 48] = true;
			this.state[xcenter + 46][ycenter + 49] = true;
			this.state[xcenter + 46][ycenter + 105] = true;
			this.state[xcenter + 46][ycenter + 116] = true;
			this.state[xcenter + 46][ycenter + 117] = true;
			this.state[xcenter + 47][ycenter + 9] = true;
			this.state[xcenter + 47][ycenter + 12] = true;
			this.state[xcenter + 47][ycenter + 13] = true;
			this.state[xcenter + 47][ycenter + 18] = true;
			this.state[xcenter + 47][ycenter + 19] = true;
			this.state[xcenter + 47][ycenter + 20] = true;
			this.state[xcenter + 47][ycenter + 24] = true;
			this.state[xcenter + 47][ycenter + 25] = true;
			this.state[xcenter + 47][ycenter + 26] = true;
			this.state[xcenter + 47][ycenter + 45] = true;
			this.state[xcenter + 47][ycenter + 49] = true;
			this.state[xcenter + 47][ycenter + 79] = true;
			this.state[xcenter + 47][ycenter + 81] = true;
			this.state[xcenter + 47][ycenter + 115] = true;
			this.state[xcenter + 47][ycenter + 116] = true;
			this.state[xcenter + 47][ycenter + 118] = true;
			this.state[xcenter + 47][ycenter + 119] = true;
			this.state[xcenter + 47][ycenter + 126] = true;
			this.state[xcenter + 47][ycenter + 127] = true;
			this.state[xcenter + 48][ycenter + 9] = true;
			this.state[xcenter + 48][ycenter + 12] = true;
			this.state[xcenter + 48][ycenter + 13] = true;
			this.state[xcenter + 48][ycenter + 19] = true;
			this.state[xcenter + 48][ycenter + 25] = true;
			this.state[xcenter + 48][ycenter + 78] = true;
			this.state[xcenter + 48][ycenter + 79] = true;
			this.state[xcenter + 48][ycenter + 81] = true;
			this.state[xcenter + 48][ycenter + 82] = true;
			this.state[xcenter + 48][ycenter + 102] = true;
			this.state[xcenter + 48][ycenter + 104] = true;
			this.state[xcenter + 48][ycenter + 106] = true;
			this.state[xcenter + 48][ycenter + 108] = true;
			this.state[xcenter + 48][ycenter + 115] = true;
			this.state[xcenter + 48][ycenter + 116] = true;
			this.state[xcenter + 48][ycenter + 118] = true;
			this.state[xcenter + 48][ycenter + 119] = true;
			this.state[xcenter + 48][ycenter + 126] = true;
			this.state[xcenter + 48][ycenter + 127] = true;
			this.state[xcenter + 48][ycenter + 128] = true;
			this.state[xcenter + 49][ycenter + 12] = true;
			this.state[xcenter + 49][ycenter + 78] = true;
			this.state[xcenter + 49][ycenter + 82] = true;
			this.state[xcenter + 49][ycenter + 101] = true;
			this.state[xcenter + 49][ycenter + 109] = true;
			this.state[xcenter + 49][ycenter + 115] = true;
			this.state[xcenter + 49][ycenter + 116] = true;
			this.state[xcenter + 49][ycenter + 118] = true;
			this.state[xcenter + 49][ycenter + 119] = true;
			this.state[xcenter + 49][ycenter + 125] = true;
			this.state[xcenter + 49][ycenter + 127] = true;
			this.state[xcenter + 49][ycenter + 128] = true;
			this.state[xcenter + 50][ycenter + 8] = true;
			this.state[xcenter + 50][ycenter + 47] = true;
			this.state[xcenter + 50][ycenter + 101] = true;
			this.state[xcenter + 50][ycenter + 109] = true;
			this.state[xcenter + 50][ycenter + 117] = true;
			this.state[xcenter + 50][ycenter + 125] = true;
			this.state[xcenter + 50][ycenter + 126] = true;
			this.state[xcenter + 50][ycenter + 127] = true;
			this.state[xcenter + 51][ycenter + 8] = true;
			this.state[xcenter + 51][ycenter + 11] = true;
			this.state[xcenter + 51][ycenter + 16] = true;
			this.state[xcenter + 51][ycenter + 17] = true;
			this.state[xcenter + 51][ycenter + 47] = true;
			this.state[xcenter + 51][ycenter + 101] = true;
			this.state[xcenter + 51][ycenter + 104] = true;
			this.state[xcenter + 51][ycenter + 106] = true;
			this.state[xcenter + 51][ycenter + 109] = true;
			this.state[xcenter + 51][ycenter + 126] = true;
			this.state[xcenter + 52][ycenter + 8] = true;
			this.state[xcenter + 52][ycenter + 10] = true;
			this.state[xcenter + 52][ycenter + 15] = true;
			this.state[xcenter + 52][ycenter + 16] = true;
			this.state[xcenter + 52][ycenter + 17] = true;
			this.state[xcenter + 52][ycenter + 46] = true;
			this.state[xcenter + 52][ycenter + 48] = true;
			this.state[xcenter + 52][ycenter + 80] = true;
			this.state[xcenter + 52][ycenter + 101] = true;
			this.state[xcenter + 52][ycenter + 102] = true;
			this.state[xcenter + 52][ycenter + 103] = true;
			this.state[xcenter + 52][ycenter + 107] = true;
			this.state[xcenter + 52][ycenter + 108] = true;
			this.state[xcenter + 52][ycenter + 109] = true;
			this.state[xcenter + 53][ycenter + 9] = true;
			this.state[xcenter + 53][ycenter + 10] = true;
			this.state[xcenter + 53][ycenter + 15] = true;
			this.state[xcenter + 53][ycenter + 16] = true;
			this.state[xcenter + 53][ycenter + 18] = true;
			this.state[xcenter + 53][ycenter + 46] = true;
			this.state[xcenter + 53][ycenter + 47] = true;
			this.state[xcenter + 53][ycenter + 48] = true;
			this.state[xcenter + 53][ycenter + 80] = true;
			this.state[xcenter + 53][ycenter + 119] = true;
			this.state[xcenter + 53][ycenter + 120] = true;
			this.state[xcenter + 54][ycenter + 8] = true;
			this.state[xcenter + 54][ycenter + 9] = true;
			this.state[xcenter + 54][ycenter + 10] = true;
			this.state[xcenter + 54][ycenter + 16] = true;
			this.state[xcenter + 54][ycenter + 17] = true;
			this.state[xcenter + 54][ycenter + 18] = true;
			this.state[xcenter + 54][ycenter + 46] = true;
			this.state[xcenter + 54][ycenter + 47] = true;
			this.state[xcenter + 54][ycenter + 48] = true;
			this.state[xcenter + 54][ycenter + 79] = true;
			this.state[xcenter + 54][ycenter + 81] = true;
			this.state[xcenter + 54][ycenter + 119] = true;
			this.state[xcenter + 54][ycenter + 120] = true;
			this.state[xcenter + 55][ycenter + 8] = true;
			this.state[xcenter + 55][ycenter + 9] = true;
			this.state[xcenter + 55][ycenter + 10] = true;
			this.state[xcenter + 55][ycenter + 17] = true;
			this.state[xcenter + 55][ycenter + 79] = true;
			this.state[xcenter + 55][ycenter + 80] = true;
			this.state[xcenter + 55][ycenter + 81] = true;
			this.state[xcenter + 55][ycenter + 110] = true;
			this.state[xcenter + 55][ycenter + 112] = true;
			this.state[xcenter + 56][ycenter + 7] = true;
			this.state[xcenter + 56][ycenter + 9] = true;
			this.state[xcenter + 56][ycenter + 10] = true;
			this.state[xcenter + 56][ycenter + 36] = true;
			this.state[xcenter + 56][ycenter + 42] = true;
			this.state[xcenter + 56][ycenter + 44] = true;
			this.state[xcenter + 56][ycenter + 50] = true;
			this.state[xcenter + 56][ycenter + 52] = true;
			this.state[xcenter + 56][ycenter + 79] = true;
			this.state[xcenter + 56][ycenter + 80] = true;
			this.state[xcenter + 56][ycenter + 81] = true;
			this.state[xcenter + 56][ycenter + 109] = true;
			this.state[xcenter + 57][ycenter + 11] = true;
			this.state[xcenter + 57][ycenter + 36] = true;
			this.state[xcenter + 57][ycenter + 37] = true;
			this.state[xcenter + 57][ycenter + 45] = true;
			this.state[xcenter + 57][ycenter + 49] = true;
			this.state[xcenter + 57][ycenter + 109] = true;
			this.state[xcenter + 58][ycenter + 7] = true;
			this.state[xcenter + 58][ycenter + 9] = true;
			this.state[xcenter + 58][ycenter + 12] = true;
			this.state[xcenter + 58][ycenter + 13] = true;
			this.state[xcenter + 58][ycenter + 35] = true;
			this.state[xcenter + 58][ycenter + 37] = true;
			this.state[xcenter + 58][ycenter + 45] = true;
			this.state[xcenter + 58][ycenter + 49] = true;
			this.state[xcenter + 58][ycenter + 75] = true;
			this.state[xcenter + 58][ycenter + 77] = true;
			this.state[xcenter + 58][ycenter + 83] = true;
			this.state[xcenter + 58][ycenter + 85] = true;
			this.state[xcenter + 58][ycenter + 91] = true;
			this.state[xcenter + 58][ycenter + 109] = true;
			this.state[xcenter + 58][ycenter + 112] = true;
			this.state[xcenter + 59][ycenter + 1] = true;
			this.state[xcenter + 59][ycenter + 3] = true;
			this.state[xcenter + 59][ycenter + 8] = true;
			this.state[xcenter + 59][ycenter + 12] = true;
			this.state[xcenter + 59][ycenter + 15] = true;
			this.state[xcenter + 59][ycenter + 17] = true;
			this.state[xcenter + 59][ycenter + 42] = true;
			this.state[xcenter + 59][ycenter + 45] = true;
			this.state[xcenter + 59][ycenter + 49] = true;
			this.state[xcenter + 59][ycenter + 52] = true;
			this.state[xcenter + 59][ycenter + 78] = true;
			this.state[xcenter + 59][ycenter + 82] = true;
			this.state[xcenter + 59][ycenter + 90] = true;
			this.state[xcenter + 59][ycenter + 91] = true;
			this.state[xcenter + 59][ycenter + 109] = true;
			this.state[xcenter + 59][ycenter + 110] = true;
			this.state[xcenter + 59][ycenter + 111] = true;
			this.state[xcenter + 60][ycenter + 4] = true;
			this.state[xcenter + 60][ycenter + 9] = true;
			this.state[xcenter + 60][ycenter + 12] = true;
			this.state[xcenter + 60][ycenter + 18] = true;
			this.state[xcenter + 60][ycenter + 43] = true;
			this.state[xcenter + 60][ycenter + 44] = true;
			this.state[xcenter + 60][ycenter + 45] = true;
			this.state[xcenter + 60][ycenter + 49] = true;
			this.state[xcenter + 60][ycenter + 50] = true;
			this.state[xcenter + 60][ycenter + 51] = true;
			this.state[xcenter + 60][ycenter + 78] = true;
			this.state[xcenter + 60][ycenter + 82] = true;
			this.state[xcenter + 60][ycenter + 90] = true;
			this.state[xcenter + 60][ycenter + 92] = true;
			this.state[xcenter + 61][ycenter + 4] = true;
			this.state[xcenter + 61][ycenter + 9] = true;
			this.state[xcenter + 61][ycenter + 10] = true;
			this.state[xcenter + 61][ycenter + 11] = true;
			this.state[xcenter + 61][ycenter + 18] = true;
			this.state[xcenter + 61][ycenter + 75] = true;
			this.state[xcenter + 61][ycenter + 78] = true;
			this.state[xcenter + 61][ycenter + 82] = true;
			this.state[xcenter + 61][ycenter + 85] = true;
			this.state[xcenter + 62][ycenter + 1] = true;
			this.state[xcenter + 62][ycenter + 4] = true;
			this.state[xcenter + 62][ycenter + 15] = true;
			this.state[xcenter + 62][ycenter + 18] = true;
			this.state[xcenter + 62][ycenter + 76] = true;
			this.state[xcenter + 62][ycenter + 77] = true;
			this.state[xcenter + 62][ycenter + 78] = true;
			this.state[xcenter + 62][ycenter + 82] = true;
			this.state[xcenter + 62][ycenter + 83] = true;
			this.state[xcenter + 62][ycenter + 84] = true;
			this.state[xcenter + 63][ycenter + 2] = true;
			this.state[xcenter + 63][ycenter + 3] = true;
			this.state[xcenter + 63][ycenter + 4] = true;
			this.state[xcenter + 63][ycenter + 16] = true;
			this.state[xcenter + 63][ycenter + 17] = true;
			this.state[xcenter + 63][ycenter + 18] = true;
			this.state[xcenter + 63][ycenter + 39] = true;
			this.state[xcenter + 63][ycenter + 41] = true;
			this.state[xcenter + 63][ycenter + 58] = true;
			this.state[xcenter + 63][ycenter + 115] = true;
			this.state[xcenter + 63][ycenter + 116] = true;
			this.state[xcenter + 63][ycenter + 117] = true;
			this.state[xcenter + 64][ycenter + 38] = true;
			this.state[xcenter + 64][ycenter + 58] = true;
			this.state[xcenter + 64][ycenter + 109] = true;
			this.state[xcenter + 64][ycenter + 110] = true;
			this.state[xcenter + 64][ycenter + 115] = true;
			this.state[xcenter + 64][ycenter + 118] = true;
			this.state[xcenter + 64][ycenter + 119] = true;
			this.state[xcenter + 64][ycenter + 123] = true;
			this.state[xcenter + 64][ycenter + 124] = true;
			this.state[xcenter + 65][ycenter + 38] = true;
			this.state[xcenter + 65][ycenter + 58] = true;
			this.state[xcenter + 65][ycenter + 69] = true;
			this.state[xcenter + 65][ycenter + 86] = true;
			this.state[xcenter + 65][ycenter + 88] = true;
			this.state[xcenter + 65][ycenter + 108] = true;
			this.state[xcenter + 65][ycenter + 109] = true;
			this.state[xcenter + 65][ycenter + 110] = true;
			this.state[xcenter + 65][ycenter + 116] = true;
			this.state[xcenter + 65][ycenter + 117] = true;
			this.state[xcenter + 65][ycenter + 118] = true;
			this.state[xcenter + 65][ycenter + 122] = true;
			this.state[xcenter + 65][ycenter + 123] = true;
			this.state[xcenter + 65][ycenter + 124] = true;
			this.state[xcenter + 66][ycenter + 38] = true;
			this.state[xcenter + 66][ycenter + 41] = true;
			this.state[xcenter + 66][ycenter + 49] = true;
			this.state[xcenter + 66][ycenter + 50] = true;
			this.state[xcenter + 66][ycenter + 69] = true;
			this.state[xcenter + 66][ycenter + 89] = true;
			this.state[xcenter + 66][ycenter + 108] = true;
			this.state[xcenter + 66][ycenter + 109] = true;
			this.state[xcenter + 66][ycenter + 111] = true;
			this.state[xcenter + 66][ycenter + 117] = true;
			this.state[xcenter + 66][ycenter + 122] = true;
			this.state[xcenter + 66][ycenter + 123] = true;
			this.state[xcenter + 66][ycenter + 125] = true;
			this.state[xcenter + 67][ycenter + 28] = true;
			this.state[xcenter + 67][ycenter + 29] = true;
			this.state[xcenter + 67][ycenter + 30] = true;
			this.state[xcenter + 67][ycenter + 38] = true;
			this.state[xcenter + 67][ycenter + 39] = true;
			this.state[xcenter + 67][ycenter + 40] = true;
			this.state[xcenter + 67][ycenter + 48] = true;
			this.state[xcenter + 67][ycenter + 49] = true;
			this.state[xcenter + 67][ycenter + 50] = true;
			this.state[xcenter + 67][ycenter + 69] = true;
			this.state[xcenter + 67][ycenter + 89] = true;
			this.state[xcenter + 67][ycenter + 109] = true;
			this.state[xcenter + 67][ycenter + 110] = true;
			this.state[xcenter + 67][ycenter + 111] = true;
			this.state[xcenter + 67][ycenter + 123] = true;
			this.state[xcenter + 67][ycenter + 124] = true;
			this.state[xcenter + 67][ycenter + 125] = true;
			this.state[xcenter + 68][ycenter + 28] = true;
			this.state[xcenter + 68][ycenter + 29] = true;
			this.state[xcenter + 68][ycenter + 30] = true;
			this.state[xcenter + 68][ycenter + 48] = true;
			this.state[xcenter + 68][ycenter + 49] = true;
			this.state[xcenter + 68][ycenter + 50] = true;
			this.state[xcenter + 68][ycenter + 56] = true;
			this.state[xcenter + 68][ycenter + 57] = true;
			this.state[xcenter + 68][ycenter + 77] = true;
			this.state[xcenter + 68][ycenter + 78] = true;
			this.state[xcenter + 68][ycenter + 86] = true;
			this.state[xcenter + 68][ycenter + 89] = true;
			this.state[xcenter + 68][ycenter + 110] = true;
			this.state[xcenter + 68][ycenter + 124] = true;
			this.state[xcenter + 69][ycenter + 27] = true;
			this.state[xcenter + 69][ycenter + 28] = true;
			this.state[xcenter + 69][ycenter + 31] = true;
			this.state[xcenter + 69][ycenter + 48] = true;
			this.state[xcenter + 69][ycenter + 49] = true;
			this.state[xcenter + 69][ycenter + 51] = true;
			this.state[xcenter + 69][ycenter + 55] = true;
			this.state[xcenter + 69][ycenter + 58] = true;
			this.state[xcenter + 69][ycenter + 77] = true;
			this.state[xcenter + 69][ycenter + 78] = true;
			this.state[xcenter + 69][ycenter + 79] = true;
			this.state[xcenter + 69][ycenter + 87] = true;
			this.state[xcenter + 69][ycenter + 88] = true;
			this.state[xcenter + 69][ycenter + 89] = true;
			this.state[xcenter + 69][ycenter + 97] = true;
			this.state[xcenter + 69][ycenter + 98] = true;
			this.state[xcenter + 69][ycenter + 99] = true;
			this.state[xcenter + 70][ycenter + 27] = true;
			this.state[xcenter + 70][ycenter + 28] = true;
			this.state[xcenter + 70][ycenter + 31] = true;
			this.state[xcenter + 70][ycenter + 49] = true;
			this.state[xcenter + 70][ycenter + 50] = true;
			this.state[xcenter + 70][ycenter + 51] = true;
			this.state[xcenter + 70][ycenter + 55] = true;
			this.state[xcenter + 70][ycenter + 57] = true;
			this.state[xcenter + 70][ycenter + 70] = true;
			this.state[xcenter + 70][ycenter + 71] = true;
			this.state[xcenter + 70][ycenter + 77] = true;
			this.state[xcenter + 70][ycenter + 78] = true;
			this.state[xcenter + 70][ycenter + 79] = true;
			this.state[xcenter + 70][ycenter + 97] = true;
			this.state[xcenter + 70][ycenter + 98] = true;
			this.state[xcenter + 70][ycenter + 99] = true;
			this.state[xcenter + 71][ycenter + 28] = true;
			this.state[xcenter + 71][ycenter + 50] = true;
			this.state[xcenter + 71][ycenter + 56] = true;
			this.state[xcenter + 71][ycenter + 69] = true;
			this.state[xcenter + 71][ycenter + 72] = true;
			this.state[xcenter + 71][ycenter + 76] = true;
			this.state[xcenter + 71][ycenter + 78] = true;
			this.state[xcenter + 71][ycenter + 79] = true;
			this.state[xcenter + 71][ycenter + 96] = true;
			this.state[xcenter + 71][ycenter + 99] = true;
			this.state[xcenter + 71][ycenter + 100] = true;
			this.state[xcenter + 72][ycenter + 32] = true;
			this.state[xcenter + 72][ycenter + 70] = true;
			this.state[xcenter + 72][ycenter + 72] = true;
			this.state[xcenter + 72][ycenter + 76] = true;
			this.state[xcenter + 72][ycenter + 77] = true;
			this.state[xcenter + 72][ycenter + 78] = true;
			this.state[xcenter + 72][ycenter + 96] = true;
			this.state[xcenter + 72][ycenter + 99] = true;
			this.state[xcenter + 72][ycenter + 100] = true;
			this.state[xcenter + 73][ycenter + 23] = true;
			this.state[xcenter + 73][ycenter + 24] = true;
			this.state[xcenter + 73][ycenter + 29] = true;
			this.state[xcenter + 73][ycenter + 32] = true;
			this.state[xcenter + 73][ycenter + 57] = true;
			this.state[xcenter + 73][ycenter + 71] = true;
			this.state[xcenter + 73][ycenter + 77] = true;
			this.state[xcenter + 73][ycenter + 99] = true;
			this.state[xcenter + 74][ycenter + 23] = true;
			this.state[xcenter + 74][ycenter + 24] = true;
			this.state[xcenter + 74][ycenter + 25] = true;
			this.state[xcenter + 74][ycenter + 30] = true;
			this.state[xcenter + 74][ycenter + 32] = true;
			this.state[xcenter + 74][ycenter + 57] = true;
			this.state[xcenter + 74][ycenter + 95] = true;
			this.state[xcenter + 75][ycenter + 22] = true;
			this.state[xcenter + 75][ycenter + 24] = true;
			this.state[xcenter + 75][ycenter + 25] = true;
			this.state[xcenter + 75][ycenter + 30] = true;
			this.state[xcenter + 75][ycenter + 31] = true;
			this.state[xcenter + 75][ycenter + 70] = true;
			this.state[xcenter + 75][ycenter + 95] = true;
			this.state[xcenter + 75][ycenter + 98] = true;
			this.state[xcenter + 75][ycenter + 103] = true;
			this.state[xcenter + 75][ycenter + 104] = true;
			this.state[xcenter + 76][ycenter + 22] = true;
			this.state[xcenter + 76][ycenter + 23] = true;
			this.state[xcenter + 76][ycenter + 24] = true;
			this.state[xcenter + 76][ycenter + 30] = true;
			this.state[xcenter + 76][ycenter + 31] = true;
			this.state[xcenter + 76][ycenter + 32] = true;
			this.state[xcenter + 76][ycenter + 70] = true;
			this.state[xcenter + 76][ycenter + 95] = true;
			this.state[xcenter + 76][ycenter + 97] = true;
			this.state[xcenter + 76][ycenter + 102] = true;
			this.state[xcenter + 76][ycenter + 103] = true;
			this.state[xcenter + 76][ycenter + 104] = true;
			this.state[xcenter + 77][ycenter + 23] = true;
			this.state[xcenter + 77][ycenter + 30] = true;
			this.state[xcenter + 77][ycenter + 31] = true;
			this.state[xcenter + 77][ycenter + 32] = true;
			this.state[xcenter + 77][ycenter + 96] = true;
			this.state[xcenter + 77][ycenter + 97] = true;
			this.state[xcenter + 77][ycenter + 102] = true;
			this.state[xcenter + 77][ycenter + 103] = true;
			this.state[xcenter + 77][ycenter + 105] = true;
			this.state[xcenter + 78][ycenter + 30] = true;
			this.state[xcenter + 78][ycenter + 31] = true;
			this.state[xcenter + 78][ycenter + 33] = true;
			this.state[xcenter + 78][ycenter + 54] = true;
			this.state[xcenter + 78][ycenter + 55] = true;
			this.state[xcenter + 78][ycenter + 59] = true;
			this.state[xcenter + 78][ycenter + 60] = true;
			this.state[xcenter + 78][ycenter + 95] = true;
			this.state[xcenter + 78][ycenter + 96] = true;
			this.state[xcenter + 78][ycenter + 97] = true;
			this.state[xcenter + 78][ycenter + 103] = true;
			this.state[xcenter + 78][ycenter + 104] = true;
			this.state[xcenter + 78][ycenter + 105] = true;
			this.state[xcenter + 79][ycenter + 29] = true;
			this.state[xcenter + 79][ycenter + 95] = true;
			this.state[xcenter + 79][ycenter + 96] = true;
			this.state[xcenter + 79][ycenter + 97] = true;
			this.state[xcenter + 79][ycenter + 104] = true;
			this.state[xcenter + 80][ycenter + 27] = true;
			this.state[xcenter + 80][ycenter + 28] = true;
			this.state[xcenter + 80][ycenter + 31] = true;
			this.state[xcenter + 80][ycenter + 33] = true;
			this.state[xcenter + 80][ycenter + 67] = true;
			this.state[xcenter + 80][ycenter + 68] = true;
			this.state[xcenter + 80][ycenter + 72] = true;
			this.state[xcenter + 80][ycenter + 73] = true;
			this.state[xcenter + 80][ycenter + 94] = true;
			this.state[xcenter + 80][ycenter + 96] = true;
			this.state[xcenter + 80][ycenter + 97] = true;
			this.state[xcenter + 81][ycenter + 23] = true;
			this.state[xcenter + 81][ycenter + 25] = true;
			this.state[xcenter + 81][ycenter + 28] = true;
			this.state[xcenter + 81][ycenter + 32] = true;
			this.state[xcenter + 81][ycenter + 37] = true;
			this.state[xcenter + 81][ycenter + 39] = true;
			this.state[xcenter + 81][ycenter + 57] = true;
			this.state[xcenter + 81][ycenter + 98] = true;
			this.state[xcenter + 82][ycenter + 22] = true;
			this.state[xcenter + 82][ycenter + 28] = true;
			this.state[xcenter + 82][ycenter + 31] = true;
			this.state[xcenter + 82][ycenter + 36] = true;
			this.state[xcenter + 82][ycenter + 51] = true;
			this.state[xcenter + 82][ycenter + 53] = true;
			this.state[xcenter + 82][ycenter + 56] = true;
			this.state[xcenter + 82][ycenter + 58] = true;
			this.state[xcenter + 82][ycenter + 94] = true;
			this.state[xcenter + 82][ycenter + 96] = true;
			this.state[xcenter + 82][ycenter + 99] = true;
			this.state[xcenter + 82][ycenter + 100] = true;
			this.state[xcenter + 83][ycenter + 22] = true;
			this.state[xcenter + 83][ycenter + 29] = true;
			this.state[xcenter + 83][ycenter + 30] = true;
			this.state[xcenter + 83][ycenter + 31] = true;
			this.state[xcenter + 83][ycenter + 36] = true;
			this.state[xcenter + 83][ycenter + 52] = true;
			this.state[xcenter + 83][ycenter + 53] = true;
			this.state[xcenter + 83][ycenter + 56] = true;
			this.state[xcenter + 83][ycenter + 58] = true;
			this.state[xcenter + 83][ycenter + 70] = true;
			this.state[xcenter + 83][ycenter + 88] = true;
			this.state[xcenter + 83][ycenter + 90] = true;
			this.state[xcenter + 83][ycenter + 95] = true;
			this.state[xcenter + 83][ycenter + 99] = true;
			this.state[xcenter + 83][ycenter + 102] = true;
			this.state[xcenter + 83][ycenter + 104] = true;
			this.state[xcenter + 84][ycenter + 22] = true;
			this.state[xcenter + 84][ycenter + 25] = true;
			this.state[xcenter + 84][ycenter + 36] = true;
			this.state[xcenter + 84][ycenter + 39] = true;
			this.state[xcenter + 84][ycenter + 52] = true;
			this.state[xcenter + 84][ycenter + 57] = true;
			this.state[xcenter + 84][ycenter + 69] = true;
			this.state[xcenter + 84][ycenter + 71] = true;
			this.state[xcenter + 84][ycenter + 74] = true;
			this.state[xcenter + 84][ycenter + 76] = true;
			this.state[xcenter + 84][ycenter + 91] = true;
			this.state[xcenter + 84][ycenter + 96] = true;
			this.state[xcenter + 84][ycenter + 99] = true;
			this.state[xcenter + 84][ycenter + 105] = true;
			this.state[xcenter + 85][ycenter + 22] = true;
			this.state[xcenter + 85][ycenter + 23] = true;
			this.state[xcenter + 85][ycenter + 24] = true;
			this.state[xcenter + 85][ycenter + 36] = true;
			this.state[xcenter + 85][ycenter + 37] = true;
			this.state[xcenter + 85][ycenter + 38] = true;
			this.state[xcenter + 85][ycenter + 57] = true;
			this.state[xcenter + 85][ycenter + 69] = true;
			this.state[xcenter + 85][ycenter + 71] = true;
			this.state[xcenter + 85][ycenter + 74] = true;
			this.state[xcenter + 85][ycenter + 75] = true;
			this.state[xcenter + 85][ycenter + 91] = true;
			this.state[xcenter + 85][ycenter + 96] = true;
			this.state[xcenter + 85][ycenter + 97] = true;
			this.state[xcenter + 85][ycenter + 98] = true;
			this.state[xcenter + 85][ycenter + 105] = true;
			this.state[xcenter + 86][ycenter + 56] = true;
			this.state[xcenter + 86][ycenter + 58] = true;
			this.state[xcenter + 86][ycenter + 70] = true;
			this.state[xcenter + 86][ycenter + 75] = true;
			this.state[xcenter + 86][ycenter + 88] = true;
			this.state[xcenter + 86][ycenter + 91] = true;
			this.state[xcenter + 86][ycenter + 102] = true;
			this.state[xcenter + 86][ycenter + 105] = true;
			this.state[xcenter + 87][ycenter + 46] = true;
			this.state[xcenter + 87][ycenter + 48] = true;
			this.state[xcenter + 87][ycenter + 55] = true;
			this.state[xcenter + 87][ycenter + 59] = true;
			this.state[xcenter + 87][ycenter + 70] = true;
			this.state[xcenter + 87][ycenter + 89] = true;
			this.state[xcenter + 87][ycenter + 90] = true;
			this.state[xcenter + 87][ycenter + 91] = true;
			this.state[xcenter + 87][ycenter + 103] = true;
			this.state[xcenter + 87][ycenter + 104] = true;
			this.state[xcenter + 87][ycenter + 105] = true;
			this.state[xcenter + 88][ycenter + 47] = true;
			this.state[xcenter + 88][ycenter + 48] = true;
			this.state[xcenter + 88][ycenter + 56] = true;
			this.state[xcenter + 88][ycenter + 58] = true;
			this.state[xcenter + 88][ycenter + 69] = true;
			this.state[xcenter + 88][ycenter + 71] = true;
			this.state[xcenter + 89][ycenter + 47] = true;
			this.state[xcenter + 89][ycenter + 57] = true;
			this.state[xcenter + 89][ycenter + 68] = true;
			this.state[xcenter + 89][ycenter + 72] = true;
			this.state[xcenter + 89][ycenter + 79] = true;
			this.state[xcenter + 89][ycenter + 81] = true;
			this.state[xcenter + 90][ycenter + 69] = true;
			this.state[xcenter + 90][ycenter + 71] = true;
			this.state[xcenter + 90][ycenter + 79] = true;
			this.state[xcenter + 90][ycenter + 80] = true;
			this.state[xcenter + 91][ycenter + 54] = true;
			this.state[xcenter + 91][ycenter + 55] = true;
			this.state[xcenter + 91][ycenter + 59] = true;
			this.state[xcenter + 91][ycenter + 60] = true;
			this.state[xcenter + 91][ycenter + 70] = true;
			this.state[xcenter + 91][ycenter + 80] = true;
			this.state[xcenter + 92][ycenter + 54] = true;
			this.state[xcenter + 92][ycenter + 55] = true;
			this.state[xcenter + 92][ycenter + 56] = true;
			this.state[xcenter + 92][ycenter + 58] = true;
			this.state[xcenter + 92][ycenter + 59] = true;
			this.state[xcenter + 92][ycenter + 60] = true;
			this.state[xcenter + 93][ycenter + 48] = true;
			this.state[xcenter + 93][ycenter + 49] = true;
			this.state[xcenter + 93][ycenter + 53] = true;
			this.state[xcenter + 93][ycenter + 55] = true;
			this.state[xcenter + 93][ycenter + 56] = true;
			this.state[xcenter + 93][ycenter + 58] = true;
			this.state[xcenter + 93][ycenter + 59] = true;
			this.state[xcenter + 93][ycenter + 61] = true;
			this.state[xcenter + 93][ycenter + 67] = true;
			this.state[xcenter + 93][ycenter + 68] = true;
			this.state[xcenter + 93][ycenter + 72] = true;
			this.state[xcenter + 93][ycenter + 73] = true;
			this.state[xcenter + 94][ycenter + 47] = true;
			this.state[xcenter + 94][ycenter + 48] = true;
			this.state[xcenter + 94][ycenter + 49] = true;
			this.state[xcenter + 94][ycenter + 53] = true;
			this.state[xcenter + 94][ycenter + 54] = true;
			this.state[xcenter + 94][ycenter + 55] = true;
			this.state[xcenter + 94][ycenter + 59] = true;
			this.state[xcenter + 94][ycenter + 60] = true;
			this.state[xcenter + 94][ycenter + 61] = true;
			this.state[xcenter + 94][ycenter + 67] = true;
			this.state[xcenter + 94][ycenter + 68] = true;
			this.state[xcenter + 94][ycenter + 69] = true;
			this.state[xcenter + 94][ycenter + 71] = true;
			this.state[xcenter + 94][ycenter + 72] = true;
			this.state[xcenter + 94][ycenter + 73] = true;
			this.state[xcenter + 95][ycenter + 46] = true;
			this.state[xcenter + 95][ycenter + 47] = true;
			this.state[xcenter + 95][ycenter + 48] = true;
			this.state[xcenter + 95][ycenter + 49] = true;
			this.state[xcenter + 95][ycenter + 50] = true;
			this.state[xcenter + 95][ycenter + 54] = true;
			this.state[xcenter + 95][ycenter + 60] = true;
			this.state[xcenter + 95][ycenter + 66] = true;
			this.state[xcenter + 95][ycenter + 68] = true;
			this.state[xcenter + 95][ycenter + 69] = true;
			this.state[xcenter + 95][ycenter + 71] = true;
			this.state[xcenter + 95][ycenter + 72] = true;
			this.state[xcenter + 95][ycenter + 74] = true;
			this.state[xcenter + 95][ycenter + 78] = true;
			this.state[xcenter + 95][ycenter + 79] = true;
			this.state[xcenter + 96][ycenter + 39] = true;
			this.state[xcenter + 96][ycenter + 40] = true;
			this.state[xcenter + 96][ycenter + 45] = true;
			this.state[xcenter + 96][ycenter + 46] = true;
			this.state[xcenter + 96][ycenter + 66] = true;
			this.state[xcenter + 96][ycenter + 67] = true;
			this.state[xcenter + 96][ycenter + 68] = true;
			this.state[xcenter + 96][ycenter + 72] = true;
			this.state[xcenter + 96][ycenter + 73] = true;
			this.state[xcenter + 96][ycenter + 74] = true;
			this.state[xcenter + 96][ycenter + 78] = true;
			this.state[xcenter + 96][ycenter + 79] = true;
			this.state[xcenter + 96][ycenter + 80] = true;
			this.state[xcenter + 97][ycenter + 39] = true;
			this.state[xcenter + 97][ycenter + 40] = true;
			this.state[xcenter + 97][ycenter + 41] = true;
			this.state[xcenter + 97][ycenter + 46] = true;
			this.state[xcenter + 97][ycenter + 47] = true;
			this.state[xcenter + 97][ycenter + 48] = true;
			this.state[xcenter + 97][ycenter + 49] = true;
			this.state[xcenter + 97][ycenter + 50] = true;
			this.state[xcenter + 97][ycenter + 67] = true;
			this.state[xcenter + 97][ycenter + 73] = true;
			this.state[xcenter + 97][ycenter + 77] = true;
			this.state[xcenter + 97][ycenter + 78] = true;
			this.state[xcenter + 97][ycenter + 79] = true;
			this.state[xcenter + 97][ycenter + 80] = true;
			this.state[xcenter + 97][ycenter + 81] = true;
			this.state[xcenter + 98][ycenter + 38] = true;
			this.state[xcenter + 98][ycenter + 40] = true;
			this.state[xcenter + 98][ycenter + 41] = true;
			this.state[xcenter + 98][ycenter + 47] = true;
			this.state[xcenter + 98][ycenter + 48] = true;
			this.state[xcenter + 98][ycenter + 49] = true;
			this.state[xcenter + 98][ycenter + 81] = true;
			this.state[xcenter + 98][ycenter + 82] = true;
			this.state[xcenter + 98][ycenter + 87] = true;
			this.state[xcenter + 98][ycenter + 88] = true;
			this.state[xcenter + 99][ycenter + 38] = true;
			this.state[xcenter + 99][ycenter + 39] = true;
			this.state[xcenter + 99][ycenter + 40] = true;
			this.state[xcenter + 99][ycenter + 48] = true;
			this.state[xcenter + 99][ycenter + 77] = true;
			this.state[xcenter + 99][ycenter + 78] = true;
			this.state[xcenter + 99][ycenter + 79] = true;
			this.state[xcenter + 99][ycenter + 80] = true;
			this.state[xcenter + 99][ycenter + 81] = true;
			this.state[xcenter + 99][ycenter + 86] = true;
			this.state[xcenter + 99][ycenter + 87] = true;
			this.state[xcenter + 99][ycenter + 88] = true;
			this.state[xcenter + 100][ycenter + 39] = true;
			this.state[xcenter + 100][ycenter + 78] = true;
			this.state[xcenter + 100][ycenter + 79] = true;
			this.state[xcenter + 100][ycenter + 80] = true;
			this.state[xcenter + 100][ycenter + 86] = true;
			this.state[xcenter + 100][ycenter + 87] = true;
			this.state[xcenter + 100][ycenter + 89] = true;
			this.state[xcenter + 101][ycenter + 45] = true;
			this.state[xcenter + 101][ycenter + 46] = true;
			this.state[xcenter + 101][ycenter + 79] = true;
			this.state[xcenter + 101][ycenter + 87] = true;
			this.state[xcenter + 101][ycenter + 88] = true;
			this.state[xcenter + 101][ycenter + 89] = true;
			this.state[xcenter + 102][ycenter + 45] = true;
			this.state[xcenter + 102][ycenter + 46] = true;
			this.state[xcenter + 102][ycenter + 88] = true;
			this.state[xcenter + 103][ycenter + 81] = true;
			this.state[xcenter + 103][ycenter + 82] = true;
			this.state[xcenter + 104][ycenter + 55] = true;
			this.state[xcenter + 104][ycenter + 57] = true;
			this.state[xcenter + 104][ycenter + 81] = true;
			this.state[xcenter + 104][ycenter + 82] = true;
			this.state[xcenter + 105][ycenter + 54] = true;
			this.state[xcenter + 106][ycenter + 54] = true;
			this.state[xcenter + 106][ycenter + 70] = true;
			this.state[xcenter + 106][ycenter + 72] = true;
			this.state[xcenter + 107][ycenter + 54] = true;
			this.state[xcenter + 107][ycenter + 57] = true;
			this.state[xcenter + 107][ycenter + 73] = true;
			this.state[xcenter + 108][ycenter + 54] = true;
			this.state[xcenter + 108][ycenter + 55] = true;
			this.state[xcenter + 108][ycenter + 56] = true;
			this.state[xcenter + 108][ycenter + 73] = true;
			this.state[xcenter + 109][ycenter + 70] = true;
			this.state[xcenter + 109][ycenter + 73] = true;
			this.state[xcenter + 110][ycenter + 49] = true;
			this.state[xcenter + 110][ycenter + 50] = true;
			this.state[xcenter + 110][ycenter + 71] = true;
			this.state[xcenter + 110][ycenter + 72] = true;
			this.state[xcenter + 110][ycenter + 73] = true;
			this.state[xcenter + 111][ycenter + 47] = true;
			this.state[xcenter + 111][ycenter + 48] = true;
			this.state[xcenter + 111][ycenter + 50] = true;
			this.state[xcenter + 112][ycenter + 46] = true;
			this.state[xcenter + 112][ycenter + 50] = true;
			this.state[xcenter + 112][ycenter + 77] = true;
			this.state[xcenter + 112][ycenter + 78] = true;
			this.state[xcenter + 113][ycenter + 40] = true;
			this.state[xcenter + 113][ycenter + 41] = true;
			this.state[xcenter + 113][ycenter + 46] = true;
			this.state[xcenter + 113][ycenter + 49] = true;
			this.state[xcenter + 113][ycenter + 50] = true;
			this.state[xcenter + 113][ycenter + 54] = true;
			this.state[xcenter + 113][ycenter + 55] = true;
			this.state[xcenter + 113][ycenter + 77] = true;
			this.state[xcenter + 113][ycenter + 79] = true;
			this.state[xcenter + 113][ycenter + 80] = true;
			this.state[xcenter + 114][ycenter + 39] = true;
			this.state[xcenter + 114][ycenter + 40] = true;
			this.state[xcenter + 114][ycenter + 41] = true;
			this.state[xcenter + 114][ycenter + 47] = true;
			this.state[xcenter + 114][ycenter + 48] = true;
			this.state[xcenter + 114][ycenter + 49] = true;
			this.state[xcenter + 114][ycenter + 53] = true;
			this.state[xcenter + 114][ycenter + 54] = true;
			this.state[xcenter + 114][ycenter + 55] = true;
			this.state[xcenter + 114][ycenter + 77] = true;
			this.state[xcenter + 114][ycenter + 81] = true;
			this.state[xcenter + 115][ycenter + 39] = true;
			this.state[xcenter + 115][ycenter + 40] = true;
			this.state[xcenter + 115][ycenter + 42] = true;
			this.state[xcenter + 115][ycenter + 48] = true;
			this.state[xcenter + 115][ycenter + 53] = true;
			this.state[xcenter + 115][ycenter + 54] = true;
			this.state[xcenter + 115][ycenter + 56] = true;
			this.state[xcenter + 115][ycenter + 72] = true;
			this.state[xcenter + 115][ycenter + 73] = true;
			this.state[xcenter + 115][ycenter + 77] = true;
			this.state[xcenter + 115][ycenter + 78] = true;
			this.state[xcenter + 115][ycenter + 81] = true;
			this.state[xcenter + 115][ycenter + 86] = true;
			this.state[xcenter + 115][ycenter + 87] = true;
			this.state[xcenter + 116][ycenter + 40] = true;
			this.state[xcenter + 116][ycenter + 41] = true;
			this.state[xcenter + 116][ycenter + 42] = true;
			this.state[xcenter + 116][ycenter + 54] = true;
			this.state[xcenter + 116][ycenter + 55] = true;
			this.state[xcenter + 116][ycenter + 56] = true;
			this.state[xcenter + 116][ycenter + 72] = true;
			this.state[xcenter + 116][ycenter + 73] = true;
			this.state[xcenter + 116][ycenter + 74] = true;
			this.state[xcenter + 116][ycenter + 78] = true;
			this.state[xcenter + 116][ycenter + 79] = true;
			this.state[xcenter + 116][ycenter + 80] = true;
			this.state[xcenter + 116][ycenter + 86] = true;
			this.state[xcenter + 116][ycenter + 87] = true;
			this.state[xcenter + 116][ycenter + 88] = true;
			this.state[xcenter + 117][ycenter + 41] = true;
			this.state[xcenter + 117][ycenter + 55] = true;
			this.state[xcenter + 117][ycenter + 71] = true;
			this.state[xcenter + 117][ycenter + 73] = true;
			this.state[xcenter + 117][ycenter + 74] = true;
			this.state[xcenter + 117][ycenter + 79] = true;
			this.state[xcenter + 117][ycenter + 85] = true;
			this.state[xcenter + 117][ycenter + 87] = true;
			this.state[xcenter + 117][ycenter + 88] = true;
			this.state[xcenter + 118][ycenter + 71] = true;
			this.state[xcenter + 118][ycenter + 72] = true;
			this.state[xcenter + 118][ycenter + 73] = true;
			this.state[xcenter + 118][ycenter + 85] = true;
			this.state[xcenter + 118][ycenter + 86] = true;
			this.state[xcenter + 118][ycenter + 87] = true;
			this.state[xcenter + 119][ycenter + 72] = true;
			this.state[xcenter + 119][ycenter + 86] = true;
		}
		// drop cordership.
		this.cordership = function() {
			var xcenter = Math.floor(Math.random() * (this.width * 0.7 + 1) + this.width * 0.15);
			var ycenter = Math.floor(Math.random() * (this.height * 0.8 + 1) + this.height * 0.1);
			this.state[xcenter + -88][ycenter + -4] = true;
			this.state[xcenter + -88][ycenter + -3] = true;
			this.state[xcenter + -87][ycenter + -4] = true;
			this.state[xcenter + -87][ycenter + -3] = true;
			this.state[xcenter + -86][ycenter + -4] = true;
			this.state[xcenter + -86][ycenter + -3] = true;
			this.state[xcenter + -85][ycenter + -4] = true;
			this.state[xcenter + -85][ycenter + -3] = true;
			this.state[xcenter + -84][ycenter + 0] = true;
			this.state[xcenter + -84][ycenter + 1] = true;
			this.state[xcenter + -84][ycenter + 12] = true;
			this.state[xcenter + -84][ycenter + 13] = true;
			this.state[xcenter + -83][ycenter + 0] = true;
			this.state[xcenter + -83][ycenter + 1] = true;
			this.state[xcenter + -83][ycenter + 12] = true;
			this.state[xcenter + -83][ycenter + 13] = true;
			this.state[xcenter + -82][ycenter + -2] = true;
			this.state[xcenter + -82][ycenter + -1] = true;
			this.state[xcenter + -82][ycenter + 12] = true;
			this.state[xcenter + -82][ycenter + 13] = true;
			this.state[xcenter + -82][ycenter + 16] = true;
			this.state[xcenter + -82][ycenter + 17] = true;
			this.state[xcenter + -81][ycenter + -2] = true;
			this.state[xcenter + -81][ycenter + -1] = true;
			this.state[xcenter + -81][ycenter + 12] = true;
			this.state[xcenter + -81][ycenter + 13] = true;
			this.state[xcenter + -81][ycenter + 16] = true;
			this.state[xcenter + -81][ycenter + 17] = true;
			this.state[xcenter + -80][ycenter + -4] = true;
			this.state[xcenter + -80][ycenter + -3] = true;
			this.state[xcenter + -80][ycenter + 4] = true;
			this.state[xcenter + -80][ycenter + 5] = true;
			this.state[xcenter + -80][ycenter + 10] = true;
			this.state[xcenter + -80][ycenter + 11] = true;
			this.state[xcenter + -79][ycenter + -4] = true;
			this.state[xcenter + -79][ycenter + -3] = true;
			this.state[xcenter + -79][ycenter + 4] = true;
			this.state[xcenter + -79][ycenter + 5] = true;
			this.state[xcenter + -79][ycenter + 10] = true;
			this.state[xcenter + -79][ycenter + 11] = true;
			this.state[xcenter + -78][ycenter + -2] = true;
			this.state[xcenter + -78][ycenter + -1] = true;
			this.state[xcenter + -78][ycenter + 4] = true;
			this.state[xcenter + -78][ycenter + 5] = true;
			this.state[xcenter + -78][ycenter + 8] = true;
			this.state[xcenter + -78][ycenter + 9] = true;
			this.state[xcenter + -78][ycenter + 12] = true;
			this.state[xcenter + -78][ycenter + 13] = true;
			this.state[xcenter + -78][ycenter + 14] = true;
			this.state[xcenter + -78][ycenter + 15] = true;
			this.state[xcenter + -77][ycenter + -2] = true;
			this.state[xcenter + -77][ycenter + -1] = true;
			this.state[xcenter + -77][ycenter + 4] = true;
			this.state[xcenter + -77][ycenter + 5] = true;
			this.state[xcenter + -77][ycenter + 8] = true;
			this.state[xcenter + -77][ycenter + 9] = true;
			this.state[xcenter + -77][ycenter + 12] = true;
			this.state[xcenter + -77][ycenter + 13] = true;
			this.state[xcenter + -77][ycenter + 14] = true;
			this.state[xcenter + -77][ycenter + 15] = true;
			this.state[xcenter + -76][ycenter + 8] = true;
			this.state[xcenter + -76][ycenter + 9] = true;
			this.state[xcenter + -76][ycenter + 12] = true;
			this.state[xcenter + -76][ycenter + 13] = true;
			this.state[xcenter + -76][ycenter + 14] = true;
			this.state[xcenter + -76][ycenter + 15] = true;
			this.state[xcenter + -75][ycenter + 8] = true;
			this.state[xcenter + -75][ycenter + 9] = true;
			this.state[xcenter + -75][ycenter + 12] = true;
			this.state[xcenter + -75][ycenter + 13] = true;
			this.state[xcenter + -75][ycenter + 14] = true;
			this.state[xcenter + -75][ycenter + 15] = true;
			this.state[xcenter + -74][ycenter + 36] = true;
			this.state[xcenter + -74][ycenter + 37] = true;
			this.state[xcenter + -74][ycenter + 38] = true;
			this.state[xcenter + -74][ycenter + 39] = true;
			this.state[xcenter + -73][ycenter + 36] = true;
			this.state[xcenter + -73][ycenter + 37] = true;
			this.state[xcenter + -73][ycenter + 38] = true;
			this.state[xcenter + -73][ycenter + 39] = true;
			this.state[xcenter + -72][ycenter + 36] = true;
			this.state[xcenter + -72][ycenter + 37] = true;
			this.state[xcenter + -72][ycenter + 38] = true;
			this.state[xcenter + -72][ycenter + 39] = true;
			this.state[xcenter + -71][ycenter + 36] = true;
			this.state[xcenter + -71][ycenter + 37] = true;
			this.state[xcenter + -71][ycenter + 38] = true;
			this.state[xcenter + -71][ycenter + 39] = true;
			this.state[xcenter + -58][ycenter + 52] = true;
			this.state[xcenter + -58][ycenter + 53] = true;
			this.state[xcenter + -58][ycenter + 54] = true;
			this.state[xcenter + -58][ycenter + 55] = true;
			this.state[xcenter + -57][ycenter + 52] = true;
			this.state[xcenter + -57][ycenter + 53] = true;
			this.state[xcenter + -57][ycenter + 54] = true;
			this.state[xcenter + -57][ycenter + 55] = true;
			this.state[xcenter + -56][ycenter + -36] = true;
			this.state[xcenter + -56][ycenter + -35] = true;
			this.state[xcenter + -56][ycenter + -34] = true;
			this.state[xcenter + -56][ycenter + -33] = true;
			this.state[xcenter + -56][ycenter + -28] = true;
			this.state[xcenter + -56][ycenter + -27] = true;
			this.state[xcenter + -56][ycenter + 52] = true;
			this.state[xcenter + -56][ycenter + 53] = true;
			this.state[xcenter + -56][ycenter + 54] = true;
			this.state[xcenter + -56][ycenter + 55] = true;
			this.state[xcenter + -55][ycenter + -36] = true;
			this.state[xcenter + -55][ycenter + -35] = true;
			this.state[xcenter + -55][ycenter + -34] = true;
			this.state[xcenter + -55][ycenter + -33] = true;
			this.state[xcenter + -55][ycenter + -28] = true;
			this.state[xcenter + -55][ycenter + -27] = true;
			this.state[xcenter + -55][ycenter + 52] = true;
			this.state[xcenter + -55][ycenter + 53] = true;
			this.state[xcenter + -55][ycenter + 54] = true;
			this.state[xcenter + -55][ycenter + 55] = true;
			this.state[xcenter + -54][ycenter + -30] = true;
			this.state[xcenter + -54][ycenter + -29] = true;
			this.state[xcenter + -54][ycenter + -26] = true;
			this.state[xcenter + -54][ycenter + -25] = true;
			this.state[xcenter + -53][ycenter + -30] = true;
			this.state[xcenter + -53][ycenter + -29] = true;
			this.state[xcenter + -53][ycenter + -26] = true;
			this.state[xcenter + -53][ycenter + -25] = true;
			this.state[xcenter + -52][ycenter + -32] = true;
			this.state[xcenter + -52][ycenter + -31] = true;
			this.state[xcenter + -51][ycenter + -32] = true;
			this.state[xcenter + -51][ycenter + -31] = true;
			this.state[xcenter + -48][ycenter + -28] = true;
			this.state[xcenter + -48][ycenter + -27] = true;
			this.state[xcenter + -48][ycenter + -26] = true;
			this.state[xcenter + -48][ycenter + -25] = true;
			this.state[xcenter + -47][ycenter + -28] = true;
			this.state[xcenter + -47][ycenter + -27] = true;
			this.state[xcenter + -47][ycenter + -26] = true;
			this.state[xcenter + -47][ycenter + -25] = true;
			this.state[xcenter + -44][ycenter + -26] = true;
			this.state[xcenter + -44][ycenter + -25] = true;
			this.state[xcenter + -44][ycenter + -24] = true;
			this.state[xcenter + -44][ycenter + -23] = true;
			this.state[xcenter + -43][ycenter + -26] = true;
			this.state[xcenter + -43][ycenter + -25] = true;
			this.state[xcenter + -43][ycenter + -24] = true;
			this.state[xcenter + -43][ycenter + -23] = true;
			this.state[xcenter + -42][ycenter + -28] = true;
			this.state[xcenter + -42][ycenter + -27] = true;
			this.state[xcenter + -42][ycenter + 68] = true;
			this.state[xcenter + -42][ycenter + 69] = true;
			this.state[xcenter + -42][ycenter + 70] = true;
			this.state[xcenter + -42][ycenter + 71] = true;
			this.state[xcenter + -41][ycenter + -28] = true;
			this.state[xcenter + -41][ycenter + -27] = true;
			this.state[xcenter + -41][ycenter + 68] = true;
			this.state[xcenter + -41][ycenter + 69] = true;
			this.state[xcenter + -41][ycenter + 70] = true;
			this.state[xcenter + -41][ycenter + 71] = true;
			this.state[xcenter + -40][ycenter + -32] = true;
			this.state[xcenter + -40][ycenter + -31] = true;
			this.state[xcenter + -40][ycenter + -30] = true;
			this.state[xcenter + -40][ycenter + -29] = true;
			this.state[xcenter + -40][ycenter + -26] = true;
			this.state[xcenter + -40][ycenter + -25] = true;
			this.state[xcenter + -40][ycenter + -24] = true;
			this.state[xcenter + -40][ycenter + -23] = true;
			this.state[xcenter + -40][ycenter + 68] = true;
			this.state[xcenter + -40][ycenter + 69] = true;
			this.state[xcenter + -40][ycenter + 70] = true;
			this.state[xcenter + -40][ycenter + 71] = true;
			this.state[xcenter + -39][ycenter + -32] = true;
			this.state[xcenter + -39][ycenter + -31] = true;
			this.state[xcenter + -39][ycenter + -30] = true;
			this.state[xcenter + -39][ycenter + -29] = true;
			this.state[xcenter + -39][ycenter + -26] = true;
			this.state[xcenter + -39][ycenter + -25] = true;
			this.state[xcenter + -39][ycenter + -24] = true;
			this.state[xcenter + -39][ycenter + -23] = true;
			this.state[xcenter + -39][ycenter + 68] = true;
			this.state[xcenter + -39][ycenter + 69] = true;
			this.state[xcenter + -39][ycenter + 70] = true;
			this.state[xcenter + -39][ycenter + 71] = true;
			this.state[xcenter + -38][ycenter + -26] = true;
			this.state[xcenter + -38][ycenter + -25] = true;
			this.state[xcenter + -38][ycenter + -24] = true;
			this.state[xcenter + -38][ycenter + -23] = true;
			this.state[xcenter + -37][ycenter + -26] = true;
			this.state[xcenter + -37][ycenter + -25] = true;
			this.state[xcenter + -37][ycenter + -24] = true;
			this.state[xcenter + -37][ycenter + -23] = true;
			this.state[xcenter + -36][ycenter + -56] = true;
			this.state[xcenter + -36][ycenter + -55] = true;
			this.state[xcenter + -36][ycenter + -30] = true;
			this.state[xcenter + -36][ycenter + -29] = true;
			this.state[xcenter + -35][ycenter + -56] = true;
			this.state[xcenter + -35][ycenter + -55] = true;
			this.state[xcenter + -35][ycenter + -30] = true;
			this.state[xcenter + -35][ycenter + -29] = true;
			this.state[xcenter + -34][ycenter + -56] = true;
			this.state[xcenter + -34][ycenter + -55] = true;
			this.state[xcenter + -34][ycenter + -10] = true;
			this.state[xcenter + -34][ycenter + -9] = true;
			this.state[xcenter + -34][ycenter + -8] = true;
			this.state[xcenter + -34][ycenter + -7] = true;
			this.state[xcenter + -34][ycenter + -6] = true;
			this.state[xcenter + -34][ycenter + -5] = true;
			this.state[xcenter + -33][ycenter + -56] = true;
			this.state[xcenter + -33][ycenter + -55] = true;
			this.state[xcenter + -33][ycenter + -10] = true;
			this.state[xcenter + -33][ycenter + -9] = true;
			this.state[xcenter + -33][ycenter + -8] = true;
			this.state[xcenter + -33][ycenter + -7] = true;
			this.state[xcenter + -33][ycenter + -6] = true;
			this.state[xcenter + -33][ycenter + -5] = true;
			this.state[xcenter + -32][ycenter + -52] = true;
			this.state[xcenter + -32][ycenter + -51] = true;
			this.state[xcenter + -32][ycenter + -40] = true;
			this.state[xcenter + -32][ycenter + -39] = true;
			this.state[xcenter + -32][ycenter + -12] = true;
			this.state[xcenter + -32][ycenter + -11] = true;
			this.state[xcenter + -31][ycenter + -52] = true;
			this.state[xcenter + -31][ycenter + -51] = true;
			this.state[xcenter + -31][ycenter + -40] = true;
			this.state[xcenter + -31][ycenter + -39] = true;
			this.state[xcenter + -31][ycenter + -12] = true;
			this.state[xcenter + -31][ycenter + -11] = true;
			this.state[xcenter + -30][ycenter + -54] = true;
			this.state[xcenter + -30][ycenter + -53] = true;
			this.state[xcenter + -30][ycenter + -40] = true;
			this.state[xcenter + -30][ycenter + -39] = true;
			this.state[xcenter + -30][ycenter + -36] = true;
			this.state[xcenter + -30][ycenter + -35] = true;
			this.state[xcenter + -30][ycenter + -14] = true;
			this.state[xcenter + -30][ycenter + -13] = true;
			this.state[xcenter + -30][ycenter + -4] = true;
			this.state[xcenter + -30][ycenter + -3] = true;
			this.state[xcenter + -30][ycenter + -2] = true;
			this.state[xcenter + -30][ycenter + -1] = true;
			this.state[xcenter + -29][ycenter + -54] = true;
			this.state[xcenter + -29][ycenter + -53] = true;
			this.state[xcenter + -29][ycenter + -40] = true;
			this.state[xcenter + -29][ycenter + -39] = true;
			this.state[xcenter + -29][ycenter + -36] = true;
			this.state[xcenter + -29][ycenter + -35] = true;
			this.state[xcenter + -29][ycenter + -14] = true;
			this.state[xcenter + -29][ycenter + -13] = true;
			this.state[xcenter + -29][ycenter + -4] = true;
			this.state[xcenter + -29][ycenter + -3] = true;
			this.state[xcenter + -29][ycenter + -2] = true;
			this.state[xcenter + -29][ycenter + -1] = true;
			this.state[xcenter + -28][ycenter + -56] = true;
			this.state[xcenter + -28][ycenter + -55] = true;
			this.state[xcenter + -28][ycenter + -48] = true;
			this.state[xcenter + -28][ycenter + -47] = true;
			this.state[xcenter + -28][ycenter + -42] = true;
			this.state[xcenter + -28][ycenter + -41] = true;
			this.state[xcenter + -28][ycenter + -16] = true;
			this.state[xcenter + -28][ycenter + -15] = true;
			this.state[xcenter + -28][ycenter + -8] = true;
			this.state[xcenter + -28][ycenter + -7] = true;
			this.state[xcenter + -27][ycenter + -56] = true;
			this.state[xcenter + -27][ycenter + -55] = true;
			this.state[xcenter + -27][ycenter + -48] = true;
			this.state[xcenter + -27][ycenter + -47] = true;
			this.state[xcenter + -27][ycenter + -42] = true;
			this.state[xcenter + -27][ycenter + -41] = true;
			this.state[xcenter + -27][ycenter + -16] = true;
			this.state[xcenter + -27][ycenter + -15] = true;
			this.state[xcenter + -27][ycenter + -8] = true;
			this.state[xcenter + -27][ycenter + -7] = true;
			this.state[xcenter + -26][ycenter + -54] = true;
			this.state[xcenter + -26][ycenter + -53] = true;
			this.state[xcenter + -26][ycenter + -48] = true;
			this.state[xcenter + -26][ycenter + -47] = true;
			this.state[xcenter + -26][ycenter + -44] = true;
			this.state[xcenter + -26][ycenter + -43] = true;
			this.state[xcenter + -26][ycenter + -40] = true;
			this.state[xcenter + -26][ycenter + -39] = true;
			this.state[xcenter + -26][ycenter + -38] = true;
			this.state[xcenter + -26][ycenter + -37] = true;
			this.state[xcenter + -26][ycenter + -16] = true;
			this.state[xcenter + -26][ycenter + -15] = true;
			this.state[xcenter + -26][ycenter + -10] = true;
			this.state[xcenter + -26][ycenter + -9] = true;
			this.state[xcenter + -26][ycenter + 0] = true;
			this.state[xcenter + -26][ycenter + 1] = true;
			this.state[xcenter + -26][ycenter + 84] = true;
			this.state[xcenter + -26][ycenter + 85] = true;
			this.state[xcenter + -26][ycenter + 86] = true;
			this.state[xcenter + -26][ycenter + 87] = true;
			this.state[xcenter + -25][ycenter + -54] = true;
			this.state[xcenter + -25][ycenter + -53] = true;
			this.state[xcenter + -25][ycenter + -48] = true;
			this.state[xcenter + -25][ycenter + -47] = true;
			this.state[xcenter + -25][ycenter + -44] = true;
			this.state[xcenter + -25][ycenter + -43] = true;
			this.state[xcenter + -25][ycenter + -40] = true;
			this.state[xcenter + -25][ycenter + -39] = true;
			this.state[xcenter + -25][ycenter + -38] = true;
			this.state[xcenter + -25][ycenter + -37] = true;
			this.state[xcenter + -25][ycenter + -16] = true;
			this.state[xcenter + -25][ycenter + -15] = true;
			this.state[xcenter + -25][ycenter + -10] = true;
			this.state[xcenter + -25][ycenter + -9] = true;
			this.state[xcenter + -25][ycenter + 0] = true;
			this.state[xcenter + -25][ycenter + 1] = true;
			this.state[xcenter + -25][ycenter + 84] = true;
			this.state[xcenter + -25][ycenter + 85] = true;
			this.state[xcenter + -25][ycenter + 86] = true;
			this.state[xcenter + -25][ycenter + 87] = true;
			this.state[xcenter + -24][ycenter + -44] = true;
			this.state[xcenter + -24][ycenter + -43] = true;
			this.state[xcenter + -24][ycenter + -40] = true;
			this.state[xcenter + -24][ycenter + -39] = true;
			this.state[xcenter + -24][ycenter + -38] = true;
			this.state[xcenter + -24][ycenter + -37] = true;
			this.state[xcenter + -24][ycenter + -16] = true;
			this.state[xcenter + -24][ycenter + -15] = true;
			this.state[xcenter + -24][ycenter + -8] = true;
			this.state[xcenter + -24][ycenter + -7] = true;
			this.state[xcenter + -24][ycenter + 0] = true;
			this.state[xcenter + -24][ycenter + 1] = true;
			this.state[xcenter + -24][ycenter + 66] = true;
			this.state[xcenter + -24][ycenter + 67] = true;
			this.state[xcenter + -24][ycenter + 84] = true;
			this.state[xcenter + -24][ycenter + 85] = true;
			this.state[xcenter + -24][ycenter + 86] = true;
			this.state[xcenter + -24][ycenter + 87] = true;
			this.state[xcenter + -23][ycenter + -44] = true;
			this.state[xcenter + -23][ycenter + -43] = true;
			this.state[xcenter + -23][ycenter + -40] = true;
			this.state[xcenter + -23][ycenter + -39] = true;
			this.state[xcenter + -23][ycenter + -38] = true;
			this.state[xcenter + -23][ycenter + -37] = true;
			this.state[xcenter + -23][ycenter + -16] = true;
			this.state[xcenter + -23][ycenter + -15] = true;
			this.state[xcenter + -23][ycenter + -8] = true;
			this.state[xcenter + -23][ycenter + -7] = true;
			this.state[xcenter + -23][ycenter + 0] = true;
			this.state[xcenter + -23][ycenter + 1] = true;
			this.state[xcenter + -23][ycenter + 66] = true;
			this.state[xcenter + -23][ycenter + 67] = true;
			this.state[xcenter + -23][ycenter + 84] = true;
			this.state[xcenter + -23][ycenter + 85] = true;
			this.state[xcenter + -23][ycenter + 86] = true;
			this.state[xcenter + -23][ycenter + 87] = true;
			this.state[xcenter + -22][ycenter + -14] = true;
			this.state[xcenter + -22][ycenter + -13] = true;
			this.state[xcenter + -22][ycenter + -12] = true;
			this.state[xcenter + -22][ycenter + -11] = true;
			this.state[xcenter + -22][ycenter + -8] = true;
			this.state[xcenter + -22][ycenter + -7] = true;
			this.state[xcenter + -22][ycenter + -4] = true;
			this.state[xcenter + -22][ycenter + -3] = true;
			this.state[xcenter + -22][ycenter + -2] = true;
			this.state[xcenter + -22][ycenter + -1] = true;
			this.state[xcenter + -22][ycenter + 0] = true;
			this.state[xcenter + -22][ycenter + 1] = true;
			this.state[xcenter + -22][ycenter + 64] = true;
			this.state[xcenter + -22][ycenter + 65] = true;
			this.state[xcenter + -22][ycenter + 68] = true;
			this.state[xcenter + -22][ycenter + 69] = true;
			this.state[xcenter + -21][ycenter + -14] = true;
			this.state[xcenter + -21][ycenter + -13] = true;
			this.state[xcenter + -21][ycenter + -12] = true;
			this.state[xcenter + -21][ycenter + -11] = true;
			this.state[xcenter + -21][ycenter + -8] = true;
			this.state[xcenter + -21][ycenter + -7] = true;
			this.state[xcenter + -21][ycenter + -4] = true;
			this.state[xcenter + -21][ycenter + -3] = true;
			this.state[xcenter + -21][ycenter + -2] = true;
			this.state[xcenter + -21][ycenter + -1] = true;
			this.state[xcenter + -21][ycenter + 0] = true;
			this.state[xcenter + -21][ycenter + 1] = true;
			this.state[xcenter + -21][ycenter + 64] = true;
			this.state[xcenter + -21][ycenter + 65] = true;
			this.state[xcenter + -21][ycenter + 68] = true;
			this.state[xcenter + -21][ycenter + 69] = true;
			this.state[xcenter + -20][ycenter + -8] = true;
			this.state[xcenter + -20][ycenter + -7] = true;
			this.state[xcenter + -19][ycenter + -8] = true;
			this.state[xcenter + -19][ycenter + -7] = true;
			this.state[xcenter + -18][ycenter + -6] = true;
			this.state[xcenter + -18][ycenter + -5] = true;
			this.state[xcenter + -18][ycenter + -4] = true;
			this.state[xcenter + -18][ycenter + -3] = true;
			this.state[xcenter + -18][ycenter + -2] = true;
			this.state[xcenter + -18][ycenter + -1] = true;
			this.state[xcenter + -18][ycenter + 0] = true;
			this.state[xcenter + -18][ycenter + 1] = true;
			this.state[xcenter + -18][ycenter + 64] = true;
			this.state[xcenter + -18][ycenter + 65] = true;
			this.state[xcenter + -18][ycenter + 70] = true;
			this.state[xcenter + -18][ycenter + 71] = true;
			this.state[xcenter + -17][ycenter + -6] = true;
			this.state[xcenter + -17][ycenter + -5] = true;
			this.state[xcenter + -17][ycenter + -4] = true;
			this.state[xcenter + -17][ycenter + -3] = true;
			this.state[xcenter + -17][ycenter + -2] = true;
			this.state[xcenter + -17][ycenter + -1] = true;
			this.state[xcenter + -17][ycenter + 0] = true;
			this.state[xcenter + -17][ycenter + 1] = true;
			this.state[xcenter + -17][ycenter + 64] = true;
			this.state[xcenter + -17][ycenter + 65] = true;
			this.state[xcenter + -17][ycenter + 70] = true;
			this.state[xcenter + -17][ycenter + 71] = true;
			this.state[xcenter + -16][ycenter + -28] = true;
			this.state[xcenter + -16][ycenter + -27] = true;
			this.state[xcenter + -16][ycenter + -26] = true;
			this.state[xcenter + -16][ycenter + -25] = true;
			this.state[xcenter + -16][ycenter + -24] = true;
			this.state[xcenter + -16][ycenter + -23] = true;
			this.state[xcenter + -16][ycenter + -2] = true;
			this.state[xcenter + -16][ycenter + -1] = true;
			this.state[xcenter + -16][ycenter + 0] = true;
			this.state[xcenter + -16][ycenter + 1] = true;
			this.state[xcenter + -16][ycenter + 68] = true;
			this.state[xcenter + -16][ycenter + 69] = true;
			this.state[xcenter + -16][ycenter + 70] = true;
			this.state[xcenter + -16][ycenter + 71] = true;
			this.state[xcenter + -15][ycenter + -28] = true;
			this.state[xcenter + -15][ycenter + -27] = true;
			this.state[xcenter + -15][ycenter + -26] = true;
			this.state[xcenter + -15][ycenter + -25] = true;
			this.state[xcenter + -15][ycenter + -24] = true;
			this.state[xcenter + -15][ycenter + -23] = true;
			this.state[xcenter + -15][ycenter + -2] = true;
			this.state[xcenter + -15][ycenter + -1] = true;
			this.state[xcenter + -15][ycenter + 0] = true;
			this.state[xcenter + -15][ycenter + 1] = true;
			this.state[xcenter + -15][ycenter + 68] = true;
			this.state[xcenter + -15][ycenter + 69] = true;
			this.state[xcenter + -15][ycenter + 70] = true;
			this.state[xcenter + -15][ycenter + 71] = true;
			this.state[xcenter + -14][ycenter + -30] = true;
			this.state[xcenter + -14][ycenter + -29] = true;
			this.state[xcenter + -14][ycenter + -22] = true;
			this.state[xcenter + -14][ycenter + -21] = true;
			this.state[xcenter + -14][ycenter + 70] = true;
			this.state[xcenter + -14][ycenter + 71] = true;
			this.state[xcenter + -13][ycenter + -30] = true;
			this.state[xcenter + -13][ycenter + -29] = true;
			this.state[xcenter + -13][ycenter + -22] = true;
			this.state[xcenter + -13][ycenter + -21] = true;
			this.state[xcenter + -13][ycenter + 70] = true;
			this.state[xcenter + -13][ycenter + 71] = true;
			this.state[xcenter + -12][ycenter + -32] = true;
			this.state[xcenter + -12][ycenter + -31] = true;
			this.state[xcenter + -12][ycenter + -22] = true;
			this.state[xcenter + -12][ycenter + -21] = true;
			this.state[xcenter + -11][ycenter + -32] = true;
			this.state[xcenter + -11][ycenter + -31] = true;
			this.state[xcenter + -11][ycenter + -22] = true;
			this.state[xcenter + -11][ycenter + -21] = true;
			this.state[xcenter + -10][ycenter + -34] = true;
			this.state[xcenter + -10][ycenter + -33] = true;
			this.state[xcenter + -10][ycenter + -26] = true;
			this.state[xcenter + -10][ycenter + -25] = true;
			this.state[xcenter + -9][ycenter + -34] = true;
			this.state[xcenter + -9][ycenter + -33] = true;
			this.state[xcenter + -9][ycenter + -26] = true;
			this.state[xcenter + -9][ycenter + -25] = true;
			this.state[xcenter + -8][ycenter + -34] = true;
			this.state[xcenter + -8][ycenter + -33] = true;
			this.state[xcenter + -8][ycenter + -28] = true;
			this.state[xcenter + -8][ycenter + -27] = true;
			this.state[xcenter + -8][ycenter + -24] = true;
			this.state[xcenter + -8][ycenter + -23] = true;
			this.state[xcenter + -8][ycenter + -22] = true;
			this.state[xcenter + -8][ycenter + -21] = true;
			this.state[xcenter + -8][ycenter + -20] = true;
			this.state[xcenter + -8][ycenter + -19] = true;
			this.state[xcenter + -7][ycenter + -34] = true;
			this.state[xcenter + -7][ycenter + -33] = true;
			this.state[xcenter + -7][ycenter + -28] = true;
			this.state[xcenter + -7][ycenter + -27] = true;
			this.state[xcenter + -7][ycenter + -24] = true;
			this.state[xcenter + -7][ycenter + -23] = true;
			this.state[xcenter + -7][ycenter + -22] = true;
			this.state[xcenter + -7][ycenter + -21] = true;
			this.state[xcenter + -7][ycenter + -20] = true;
			this.state[xcenter + -7][ycenter + -19] = true;
			this.state[xcenter + -6][ycenter + -34] = true;
			this.state[xcenter + -6][ycenter + -33] = true;
			this.state[xcenter + -6][ycenter + -18] = true;
			this.state[xcenter + -6][ycenter + -17] = true;
			this.state[xcenter + -5][ycenter + -34] = true;
			this.state[xcenter + -5][ycenter + -33] = true;
			this.state[xcenter + -5][ycenter + -18] = true;
			this.state[xcenter + -5][ycenter + -17] = true;
			this.state[xcenter + -4][ycenter + -88] = true;
			this.state[xcenter + -4][ycenter + -87] = true;
			this.state[xcenter + -4][ycenter + -86] = true;
			this.state[xcenter + -4][ycenter + -85] = true;
			this.state[xcenter + -4][ycenter + -80] = true;
			this.state[xcenter + -4][ycenter + -79] = true;
			this.state[xcenter + -4][ycenter + -30] = true;
			this.state[xcenter + -4][ycenter + -29] = true;
			this.state[xcenter + -4][ycenter + -22] = true;
			this.state[xcenter + -4][ycenter + -21] = true;
			this.state[xcenter + -4][ycenter + -18] = true;
			this.state[xcenter + -4][ycenter + -17] = true;
			this.state[xcenter + -4][ycenter + 64] = true;
			this.state[xcenter + -4][ycenter + 65] = true;
			this.state[xcenter + -4][ycenter + 66] = true;
			this.state[xcenter + -4][ycenter + 67] = true;
			this.state[xcenter + -3][ycenter + -88] = true;
			this.state[xcenter + -3][ycenter + -87] = true;
			this.state[xcenter + -3][ycenter + -86] = true;
			this.state[xcenter + -3][ycenter + -85] = true;
			this.state[xcenter + -3][ycenter + -80] = true;
			this.state[xcenter + -3][ycenter + -79] = true;
			this.state[xcenter + -3][ycenter + -30] = true;
			this.state[xcenter + -3][ycenter + -29] = true;
			this.state[xcenter + -3][ycenter + -22] = true;
			this.state[xcenter + -3][ycenter + -21] = true;
			this.state[xcenter + -3][ycenter + -18] = true;
			this.state[xcenter + -3][ycenter + -17] = true;
			this.state[xcenter + -3][ycenter + 64] = true;
			this.state[xcenter + -3][ycenter + 65] = true;
			this.state[xcenter + -3][ycenter + 66] = true;
			this.state[xcenter + -3][ycenter + 67] = true;
			this.state[xcenter + -2][ycenter + -82] = true;
			this.state[xcenter + -2][ycenter + -81] = true;
			this.state[xcenter + -2][ycenter + -78] = true;
			this.state[xcenter + -2][ycenter + -77] = true;
			this.state[xcenter + -2][ycenter + -30] = true;
			this.state[xcenter + -2][ycenter + -29] = true;
			this.state[xcenter + -2][ycenter + -22] = true;
			this.state[xcenter + -2][ycenter + -21] = true;
			this.state[xcenter + -2][ycenter + -18] = true;
			this.state[xcenter + -2][ycenter + -17] = true;
			this.state[xcenter + -2][ycenter + -16] = true;
			this.state[xcenter + -2][ycenter + -15] = true;
			this.state[xcenter + -2][ycenter + 68] = true;
			this.state[xcenter + -2][ycenter + 69] = true;
			this.state[xcenter + -1][ycenter + -82] = true;
			this.state[xcenter + -1][ycenter + -81] = true;
			this.state[xcenter + -1][ycenter + -78] = true;
			this.state[xcenter + -1][ycenter + -77] = true;
			this.state[xcenter + -1][ycenter + -30] = true;
			this.state[xcenter + -1][ycenter + -29] = true;
			this.state[xcenter + -1][ycenter + -22] = true;
			this.state[xcenter + -1][ycenter + -21] = true;
			this.state[xcenter + -1][ycenter + -18] = true;
			this.state[xcenter + -1][ycenter + -17] = true;
			this.state[xcenter + -1][ycenter + -16] = true;
			this.state[xcenter + -1][ycenter + -15] = true;
			this.state[xcenter + -1][ycenter + 68] = true;
			this.state[xcenter + -1][ycenter + 69] = true;
			this.state[xcenter + 0][ycenter + -84] = true;
			this.state[xcenter + 0][ycenter + -83] = true;
			this.state[xcenter + 0][ycenter + -26] = true;
			this.state[xcenter + 0][ycenter + -25] = true;
			this.state[xcenter + 0][ycenter + -24] = true;
			this.state[xcenter + 0][ycenter + -23] = true;
			this.state[xcenter + 0][ycenter + -22] = true;
			this.state[xcenter + 0][ycenter + -21] = true;
			this.state[xcenter + 0][ycenter + -18] = true;
			this.state[xcenter + 0][ycenter + -17] = true;
			this.state[xcenter + 0][ycenter + -16] = true;
			this.state[xcenter + 0][ycenter + -15] = true;
			this.state[xcenter + 0][ycenter + 64] = true;
			this.state[xcenter + 0][ycenter + 65] = true;
			this.state[xcenter + 0][ycenter + 66] = true;
			this.state[xcenter + 0][ycenter + 67] = true;
			this.state[xcenter + 1][ycenter + -84] = true;
			this.state[xcenter + 1][ycenter + -83] = true;
			this.state[xcenter + 1][ycenter + -26] = true;
			this.state[xcenter + 1][ycenter + -25] = true;
			this.state[xcenter + 1][ycenter + -24] = true;
			this.state[xcenter + 1][ycenter + -23] = true;
			this.state[xcenter + 1][ycenter + -22] = true;
			this.state[xcenter + 1][ycenter + -21] = true;
			this.state[xcenter + 1][ycenter + -18] = true;
			this.state[xcenter + 1][ycenter + -17] = true;
			this.state[xcenter + 1][ycenter + -16] = true;
			this.state[xcenter + 1][ycenter + -15] = true;
			this.state[xcenter + 1][ycenter + 64] = true;
			this.state[xcenter + 1][ycenter + 65] = true;
			this.state[xcenter + 1][ycenter + 66] = true;
			this.state[xcenter + 1][ycenter + 67] = true;
			this.state[xcenter + 4][ycenter + -80] = true;
			this.state[xcenter + 4][ycenter + -79] = true;
			this.state[xcenter + 4][ycenter + -78] = true;
			this.state[xcenter + 4][ycenter + -77] = true;
			this.state[xcenter + 5][ycenter + -80] = true;
			this.state[xcenter + 5][ycenter + -79] = true;
			this.state[xcenter + 5][ycenter + -78] = true;
			this.state[xcenter + 5][ycenter + -77] = true;
			this.state[xcenter + 8][ycenter + -78] = true;
			this.state[xcenter + 8][ycenter + -77] = true;
			this.state[xcenter + 8][ycenter + -76] = true;
			this.state[xcenter + 8][ycenter + -75] = true;
			this.state[xcenter + 9][ycenter + -78] = true;
			this.state[xcenter + 9][ycenter + -77] = true;
			this.state[xcenter + 9][ycenter + -76] = true;
			this.state[xcenter + 9][ycenter + -75] = true;
			this.state[xcenter + 10][ycenter + -80] = true;
			this.state[xcenter + 10][ycenter + -79] = true;
			this.state[xcenter + 11][ycenter + -80] = true;
			this.state[xcenter + 11][ycenter + -79] = true;
			this.state[xcenter + 12][ycenter + -84] = true;
			this.state[xcenter + 12][ycenter + -83] = true;
			this.state[xcenter + 12][ycenter + -82] = true;
			this.state[xcenter + 12][ycenter + -81] = true;
			this.state[xcenter + 12][ycenter + -78] = true;
			this.state[xcenter + 12][ycenter + -77] = true;
			this.state[xcenter + 12][ycenter + -76] = true;
			this.state[xcenter + 12][ycenter + -75] = true;
			this.state[xcenter + 12][ycenter + 30] = true;
			this.state[xcenter + 12][ycenter + 31] = true;
			this.state[xcenter + 12][ycenter + 34] = true;
			this.state[xcenter + 12][ycenter + 35] = true;
			this.state[xcenter + 12][ycenter + 48] = true;
			this.state[xcenter + 12][ycenter + 49] = true;
			this.state[xcenter + 12][ycenter + 52] = true;
			this.state[xcenter + 12][ycenter + 53] = true;
			this.state[xcenter + 13][ycenter + -84] = true;
			this.state[xcenter + 13][ycenter + -83] = true;
			this.state[xcenter + 13][ycenter + -82] = true;
			this.state[xcenter + 13][ycenter + -81] = true;
			this.state[xcenter + 13][ycenter + -78] = true;
			this.state[xcenter + 13][ycenter + -77] = true;
			this.state[xcenter + 13][ycenter + -76] = true;
			this.state[xcenter + 13][ycenter + -75] = true;
			this.state[xcenter + 13][ycenter + 30] = true;
			this.state[xcenter + 13][ycenter + 31] = true;
			this.state[xcenter + 13][ycenter + 34] = true;
			this.state[xcenter + 13][ycenter + 35] = true;
			this.state[xcenter + 13][ycenter + 48] = true;
			this.state[xcenter + 13][ycenter + 49] = true;
			this.state[xcenter + 13][ycenter + 52] = true;
			this.state[xcenter + 13][ycenter + 53] = true;
			this.state[xcenter + 14][ycenter + -78] = true;
			this.state[xcenter + 14][ycenter + -77] = true;
			this.state[xcenter + 14][ycenter + -76] = true;
			this.state[xcenter + 14][ycenter + -75] = true;
			this.state[xcenter + 14][ycenter + 28] = true;
			this.state[xcenter + 14][ycenter + 29] = true;
			this.state[xcenter + 14][ycenter + 48] = true;
			this.state[xcenter + 14][ycenter + 49] = true;
			this.state[xcenter + 14][ycenter + 52] = true;
			this.state[xcenter + 14][ycenter + 53] = true;
			this.state[xcenter + 15][ycenter + -78] = true;
			this.state[xcenter + 15][ycenter + -77] = true;
			this.state[xcenter + 15][ycenter + -76] = true;
			this.state[xcenter + 15][ycenter + -75] = true;
			this.state[xcenter + 15][ycenter + 28] = true;
			this.state[xcenter + 15][ycenter + 29] = true;
			this.state[xcenter + 15][ycenter + 48] = true;
			this.state[xcenter + 15][ycenter + 49] = true;
			this.state[xcenter + 15][ycenter + 52] = true;
			this.state[xcenter + 15][ycenter + 53] = true;
			this.state[xcenter + 16][ycenter + -82] = true;
			this.state[xcenter + 16][ycenter + -81] = true;
			this.state[xcenter + 16][ycenter + 30] = true;
			this.state[xcenter + 16][ycenter + 31] = true;
			this.state[xcenter + 16][ycenter + 36] = true;
			this.state[xcenter + 16][ycenter + 37] = true;
			this.state[xcenter + 16][ycenter + 50] = true;
			this.state[xcenter + 16][ycenter + 51] = true;
			this.state[xcenter + 17][ycenter + -82] = true;
			this.state[xcenter + 17][ycenter + -81] = true;
			this.state[xcenter + 17][ycenter + 30] = true;
			this.state[xcenter + 17][ycenter + 31] = true;
			this.state[xcenter + 17][ycenter + 36] = true;
			this.state[xcenter + 17][ycenter + 37] = true;
			this.state[xcenter + 17][ycenter + 50] = true;
			this.state[xcenter + 17][ycenter + 51] = true;
			this.state[xcenter + 18][ycenter + 34] = true;
			this.state[xcenter + 18][ycenter + 35] = true;
			this.state[xcenter + 18][ycenter + 36] = true;
			this.state[xcenter + 18][ycenter + 37] = true;
			this.state[xcenter + 18][ycenter + 38] = true;
			this.state[xcenter + 18][ycenter + 39] = true;
			this.state[xcenter + 19][ycenter + 34] = true;
			this.state[xcenter + 19][ycenter + 35] = true;
			this.state[xcenter + 19][ycenter + 36] = true;
			this.state[xcenter + 19][ycenter + 37] = true;
			this.state[xcenter + 19][ycenter + 38] = true;
			this.state[xcenter + 19][ycenter + 39] = true;
			this.state[xcenter + 28][ycenter + 14] = true;
			this.state[xcenter + 28][ycenter + 15] = true;
			this.state[xcenter + 29][ycenter + 14] = true;
			this.state[xcenter + 29][ycenter + 15] = true;
			this.state[xcenter + 30][ycenter + 12] = true;
			this.state[xcenter + 30][ycenter + 13] = true;
			this.state[xcenter + 30][ycenter + 16] = true;
			this.state[xcenter + 30][ycenter + 17] = true;
			this.state[xcenter + 31][ycenter + 12] = true;
			this.state[xcenter + 31][ycenter + 13] = true;
			this.state[xcenter + 31][ycenter + 16] = true;
			this.state[xcenter + 31][ycenter + 17] = true;
			this.state[xcenter + 34][ycenter + 12] = true;
			this.state[xcenter + 34][ycenter + 13] = true;
			this.state[xcenter + 34][ycenter + 18] = true;
			this.state[xcenter + 34][ycenter + 19] = true;
			this.state[xcenter + 35][ycenter + 12] = true;
			this.state[xcenter + 35][ycenter + 13] = true;
			this.state[xcenter + 35][ycenter + 18] = true;
			this.state[xcenter + 35][ycenter + 19] = true;
			this.state[xcenter + 36][ycenter + -74] = true;
			this.state[xcenter + 36][ycenter + -73] = true;
			this.state[xcenter + 36][ycenter + -72] = true;
			this.state[xcenter + 36][ycenter + -71] = true;
			this.state[xcenter + 36][ycenter + 16] = true;
			this.state[xcenter + 36][ycenter + 17] = true;
			this.state[xcenter + 36][ycenter + 18] = true;
			this.state[xcenter + 36][ycenter + 19] = true;
			this.state[xcenter + 37][ycenter + -74] = true;
			this.state[xcenter + 37][ycenter + -73] = true;
			this.state[xcenter + 37][ycenter + -72] = true;
			this.state[xcenter + 37][ycenter + -71] = true;
			this.state[xcenter + 37][ycenter + 16] = true;
			this.state[xcenter + 37][ycenter + 17] = true;
			this.state[xcenter + 37][ycenter + 18] = true;
			this.state[xcenter + 37][ycenter + 19] = true;
			this.state[xcenter + 38][ycenter + -74] = true;
			this.state[xcenter + 38][ycenter + -73] = true;
			this.state[xcenter + 38][ycenter + -72] = true;
			this.state[xcenter + 38][ycenter + -71] = true;
			this.state[xcenter + 38][ycenter + 18] = true;
			this.state[xcenter + 38][ycenter + 19] = true;
			this.state[xcenter + 39][ycenter + -74] = true;
			this.state[xcenter + 39][ycenter + -73] = true;
			this.state[xcenter + 39][ycenter + -72] = true;
			this.state[xcenter + 39][ycenter + -71] = true;
			this.state[xcenter + 39][ycenter + 18] = true;
			this.state[xcenter + 39][ycenter + 19] = true;
			this.state[xcenter + 48][ycenter + 12] = true;
			this.state[xcenter + 48][ycenter + 13] = true;
			this.state[xcenter + 48][ycenter + 14] = true;
			this.state[xcenter + 48][ycenter + 15] = true;
			this.state[xcenter + 49][ycenter + 12] = true;
			this.state[xcenter + 49][ycenter + 13] = true;
			this.state[xcenter + 49][ycenter + 14] = true;
			this.state[xcenter + 49][ycenter + 15] = true;
			this.state[xcenter + 50][ycenter + 16] = true;
			this.state[xcenter + 50][ycenter + 17] = true;
			this.state[xcenter + 51][ycenter + 16] = true;
			this.state[xcenter + 51][ycenter + 17] = true;
			this.state[xcenter + 52][ycenter + -58] = true;
			this.state[xcenter + 52][ycenter + -57] = true;
			this.state[xcenter + 52][ycenter + -56] = true;
			this.state[xcenter + 52][ycenter + -55] = true;
			this.state[xcenter + 52][ycenter + 12] = true;
			this.state[xcenter + 52][ycenter + 13] = true;
			this.state[xcenter + 52][ycenter + 14] = true;
			this.state[xcenter + 52][ycenter + 15] = true;
			this.state[xcenter + 53][ycenter + -58] = true;
			this.state[xcenter + 53][ycenter + -57] = true;
			this.state[xcenter + 53][ycenter + -56] = true;
			this.state[xcenter + 53][ycenter + -55] = true;
			this.state[xcenter + 53][ycenter + 12] = true;
			this.state[xcenter + 53][ycenter + 13] = true;
			this.state[xcenter + 53][ycenter + 14] = true;
			this.state[xcenter + 53][ycenter + 15] = true;
			this.state[xcenter + 54][ycenter + -58] = true;
			this.state[xcenter + 54][ycenter + -57] = true;
			this.state[xcenter + 54][ycenter + -56] = true;
			this.state[xcenter + 54][ycenter + -55] = true;
			this.state[xcenter + 55][ycenter + -58] = true;
			this.state[xcenter + 55][ycenter + -57] = true;
			this.state[xcenter + 55][ycenter + -56] = true;
			this.state[xcenter + 55][ycenter + -55] = true;
			this.state[xcenter + 64][ycenter + -22] = true;
			this.state[xcenter + 64][ycenter + -21] = true;
			this.state[xcenter + 64][ycenter + -18] = true;
			this.state[xcenter + 64][ycenter + -17] = true;
			this.state[xcenter + 64][ycenter + -4] = true;
			this.state[xcenter + 64][ycenter + -3] = true;
			this.state[xcenter + 64][ycenter + 0] = true;
			this.state[xcenter + 64][ycenter + 1] = true;
			this.state[xcenter + 65][ycenter + -22] = true;
			this.state[xcenter + 65][ycenter + -21] = true;
			this.state[xcenter + 65][ycenter + -18] = true;
			this.state[xcenter + 65][ycenter + -17] = true;
			this.state[xcenter + 65][ycenter + -4] = true;
			this.state[xcenter + 65][ycenter + -3] = true;
			this.state[xcenter + 65][ycenter + 0] = true;
			this.state[xcenter + 65][ycenter + 1] = true;
			this.state[xcenter + 66][ycenter + -24] = true;
			this.state[xcenter + 66][ycenter + -23] = true;
			this.state[xcenter + 66][ycenter + -4] = true;
			this.state[xcenter + 66][ycenter + -3] = true;
			this.state[xcenter + 66][ycenter + 0] = true;
			this.state[xcenter + 66][ycenter + 1] = true;
			this.state[xcenter + 67][ycenter + -24] = true;
			this.state[xcenter + 67][ycenter + -23] = true;
			this.state[xcenter + 67][ycenter + -4] = true;
			this.state[xcenter + 67][ycenter + -3] = true;
			this.state[xcenter + 67][ycenter + 0] = true;
			this.state[xcenter + 67][ycenter + 1] = true;
			this.state[xcenter + 68][ycenter + -42] = true;
			this.state[xcenter + 68][ycenter + -41] = true;
			this.state[xcenter + 68][ycenter + -40] = true;
			this.state[xcenter + 68][ycenter + -39] = true;
			this.state[xcenter + 68][ycenter + -22] = true;
			this.state[xcenter + 68][ycenter + -21] = true;
			this.state[xcenter + 68][ycenter + -16] = true;
			this.state[xcenter + 68][ycenter + -15] = true;
			this.state[xcenter + 68][ycenter + -2] = true;
			this.state[xcenter + 68][ycenter + -1] = true;
			this.state[xcenter + 69][ycenter + -42] = true;
			this.state[xcenter + 69][ycenter + -41] = true;
			this.state[xcenter + 69][ycenter + -40] = true;
			this.state[xcenter + 69][ycenter + -39] = true;
			this.state[xcenter + 69][ycenter + -22] = true;
			this.state[xcenter + 69][ycenter + -21] = true;
			this.state[xcenter + 69][ycenter + -16] = true;
			this.state[xcenter + 69][ycenter + -15] = true;
			this.state[xcenter + 69][ycenter + -2] = true;
			this.state[xcenter + 69][ycenter + -1] = true;
			this.state[xcenter + 70][ycenter + -42] = true;
			this.state[xcenter + 70][ycenter + -41] = true;
			this.state[xcenter + 70][ycenter + -40] = true;
			this.state[xcenter + 70][ycenter + -39] = true;
			this.state[xcenter + 70][ycenter + -18] = true;
			this.state[xcenter + 70][ycenter + -17] = true;
			this.state[xcenter + 70][ycenter + -16] = true;
			this.state[xcenter + 70][ycenter + -15] = true;
			this.state[xcenter + 70][ycenter + -14] = true;
			this.state[xcenter + 70][ycenter + -13] = true;
			this.state[xcenter + 71][ycenter + -42] = true;
			this.state[xcenter + 71][ycenter + -41] = true;
			this.state[xcenter + 71][ycenter + -40] = true;
			this.state[xcenter + 71][ycenter + -39] = true;
			this.state[xcenter + 71][ycenter + -18] = true;
			this.state[xcenter + 71][ycenter + -17] = true;
			this.state[xcenter + 71][ycenter + -16] = true;
			this.state[xcenter + 71][ycenter + -15] = true;
			this.state[xcenter + 71][ycenter + -14] = true;
			this.state[xcenter + 71][ycenter + -13] = true;
			this.state[xcenter + 84][ycenter + -26] = true;
			this.state[xcenter + 84][ycenter + -25] = true;
			this.state[xcenter + 84][ycenter + -24] = true;
			this.state[xcenter + 84][ycenter + -23] = true;
			this.state[xcenter + 85][ycenter + -26] = true;
			this.state[xcenter + 85][ycenter + -25] = true;
			this.state[xcenter + 85][ycenter + -24] = true;
			this.state[xcenter + 85][ycenter + -23] = true;
			this.state[xcenter + 86][ycenter + -26] = true;
			this.state[xcenter + 86][ycenter + -25] = true;
			this.state[xcenter + 86][ycenter + -24] = true;
			this.state[xcenter + 86][ycenter + -23] = true;
			this.state[xcenter + 87][ycenter + -26] = true;
			this.state[xcenter + 87][ycenter + -25] = true;
			this.state[xcenter + 87][ycenter + -24] = true;
			this.state[xcenter + 87][ycenter + -23] = true;
	};
}
// x,v paper points
function PondPoint (x, v) {
	this.x = x;
	this.v = v;
}

// requires paper.js
function Pond ($ctnr) {
	var self = this; // <3 JS!
	this.paperScope = new paper.PaperScope(); // so we can reuse this code
	this.paperScope.activate();
	this.INSET = 50;
	this.E = .2;
	this.$ctnr = $ctnr;
	this.width = Math.round($ctnr.width());
	this.height = Math.round($ctnr.height());
	this.canvas = document.createElement("canvas");
	$(this.canvas).css("background-color", "#006699"); // <-- key aspect of the pond-ness
	this.canvas.setAttribute("id", "pond");
	this.canvas.setAttribute("width", this.width);
	this.canvas.setAttribute("height", this.height);
	this.$ctnr.append($(this.canvas));
	this.context = this.canvas.getContext("2d");

	paper.setup("pond");

	this.ripples = []; // store paths as arrays of PondPoints, to be stepped & rendered
	this.pathsGroup = new paper.Group(); // store current rendered paths, for ease of removal

	// the shore is where the waves start coming back to the center, used in this.step
	var spoint = new paper.Point(this.INSET,this.INSET);
	var ssize = new paper.Size(self.width - this.INSET * 2, self.height - this.INSET * 2);;
	this.shore = new paper.Path.Rectangle(spoint, ssize);
	this.shore.visible = false;

	this.render = function() {
		this.pathsGroup.removeChildren();

		for (var i = 0; i < this.ripples.length; i++) {
			var path = new paper.Path();
			path.strokeColor = "white";
			for (var p = 0; p < this.ripples[i].length; p++) {
				path.add(this.ripples[i][p].x);
			}
			path.closed = true;
			path.strokeColor = "white";
			path.fillColor = "#dee";
			path.opacity = 0.2;
			path.smooth();
			this.pathsGroup.addChild(path);
		}
	};

	this.step = function(dt) {
		for (var i = 0; i < this.ripples.length; i++) {
			var pt;
			for (var p = 0; p < this.ripples[i].length; p++) {
				pt = this.ripples[i][p];
				// accelerate away from shore edges
				if (!this.shore.contains(pt.x)) {
					if (pt.x.y <= this.INSET) {
						pt.v.y += this.E;
					} else if (pt.x.y >= this.height - this.INSET * 2) {
						pt.v.y -= this.E;
					}
					if (pt.x.x <= this.INSET) {
						pt.v.x += this.E;
					} else if (pt.x.x >= this.width - this.INSET * 2) {
						pt.v.x -= this.E;
					}
					pt.v = pt.v.multiply(0.99);
				}
				pt.x = pt.x.add(pt.v.multiply(dt / 1000)); //ms
			}
		}
	};

	this.update = function() {
		if (self.ripples.length > 0) {
			var dt = Date.now() - self.time;
			self.step(dt);
			self.render();
			self.time = Date.now();
		}
	};

	this.time = Date.now();
	paper.view.on('frame', this.update);
	// setInterval(this.update, 100);	

	this.tool = new paper.Tool();
	this.tool.onMouseDown = function(event) {
		self.time = Date.now();
		self.ripples.push([]);
		var pt;
		// surround event.point with path points
		for (var theta = 0; theta < 2* Math.PI; theta += Math.PI/64) {
			pt = new paper.Point(Math.cos(theta), Math.sin(theta));
			var ppt = new PondPoint(event.point.add(pt), pt.multiply(30));
			self.ripples[self.ripples.length - 1].push(ppt);
		};
	}
}
// Lotke-Volterra predator-prey system ODE:
//   dx/dt = ax - bxy = f(x,y)
//   dy/dt = dxy - cy = g(x,y)
//   @Params: a,b,c,d,xo,yo > 0, all reals
function Model (a,b,c,d,xo,yo) {
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

/**
 * Converts an HSL color value to RGB. Conversion formula
 * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
 * Assumes h, s, and l are contained in the set [0, 1] and
 * returns r, g, and b in the set [0, 255].
 *
 * @param   Number  h       The hue
 * @param   Number  s       The saturation
 * @param   Number  l       The lightness
 * @return  Array           The RGB representation
 */
function hslToRgb(h, s, l){
    var r, g, b;

    if (s == 0) {
        r = g = b = l; // achromatic
    } else {
        var hue2rgb = function hue2rgb(p, q, t){
            if(t < 0) t += 1;
            if(t > 1) t -= 1;
            if(t < 1/6) return p + (q - p) * 6 * t;
            if(t < 1/2) return q;
            if(t < 2/3) return p + (q - p) * (2/3 - t) * 6;
            return p;
        }

        var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
        var p = 2 * l - q;
        r = hue2rgb(p, q, h + 1/3);
        g = hue2rgb(p, q, h);
        b = hue2rgb(p, q, h - 1/3);
    }

    var rs = Math.round(r * 255).toString(16);
    var gs = Math.round(g * 255).toString(16);
    var bs = Math.round(b * 255).toString(16);

    rs = (rs.length > 1) ? rs : "0" + rs;
    gs = (gs.length > 1) ? gs : "0" + gs;
    bs = (bs.length > 1) ? bs : "0" + bs;

    return "#" + rs + gs + bs;
};

var TID = undefined;
var animate = function() {
	if (TID) clearInterval(TID);

	var pp = new Model (Math.random(),Math.random(),Math.random(),Math.random(),Math.random(),Math.random());
	var $bg = $("body");

	TID = setInterval(function () {
		var c = hslToRgb(pp.x, pp.y, 0.75);
		// console.log("Color: " + c.length);
		$bg.css("background-color", c);
		movept(pp.x,pp.y);
		pp.step(.001);
	}, 1);
};
function RuleOneTen ($container) { // takes a JQ container div
	this.generation = 0;
	this.state = [];
	this.width = Math.round($container.width());
	this.height = Math.round($container.height());
	this.loopid = null;
	this.$container = $container;

	// init empty generation
	for (var iii = 0; iii < this.width; iii++) {
		this.state.push(!Math.round(Math.random()));
	};

	// init canvas
	this.canvas = document.createElement("canvas");
	this.canvas.setAttribute("id", "conway");
	this.canvas.setAttribute("width", this.width);
	this.canvas.setAttribute("height", this.height);
	this.$container.append($(this.canvas));
	this.context = this.canvas.getContext("2d");
	this.context.imageSmoothingEnabled= false;

	// survive to next g?
	// don't run edge cases here! i.e. col 0, col (width)
	this.lives = function(i) {
		var p = this.state[i - 1] || false;
		var q = this.state[i] || false;
		var r = this.state[i + 1] || false;
		return (q && !p) || ((q || r) && !(q && r)); // rule 110
	};

	// compute the next generation
	this.step = function() {
		var temp = [];
		// not edge cases
		for (var i = 0; i < this.width; i++) {
			temp.push(this.lives(i));
		}
		this.state = temp;
		this.generation++;
	};

	this.drawRow = function() {
		var r = this.generation % this.height;
		for (var i = 0; i < this.width; i++) {
			if (this.state[i]) {
				this.context.fillStyle = "#fff";
			} else {
				this.context.fillStyle = "#000";
			}
			this.context.fillRect(i,r,1,1);
		}
	};

	var self = this; // love JS!
	this.destroy = function() {
		clearInterval(self.loopid);
	};

	this.loop = function() {
		this.loopid = setInterval(function() {
			self.step();
			self.drawRow();
		}, 1);
	};

	this.reSeed = function() {
		clearInterval(this.loopid);
		this.context.clearRect(0,0,this.width,this.height);
		this.state = [];
		for (var iii = 0; iii < this.width; iii++) {
			this.state.push(!Math.round(Math.random()));
		}
		this.generation = 0;
		this.loop();
	};
};
// <script src="https://code.jquery.com/jquery-2.1.3.min.js"></script> !!
// <script src="https://code.jquery.com/ui/1.11.1/jquery-ui.min.js"></script> !!
// style it yrself!
function Window (title, content, close, width, height) {
	this.title = title;
	this.content = content;
	this.id = Date.now();
	this.width = width;
	this.height = height;
	this.close = close;
	this.onClose = null; // window close handler 4 cleanup

	this.html = "<div class='dragme' data-window-id='" + this.id + "'>" +
					"<div class='handle'>" + this.title  + (close ? "<p class='x'>&times;</p>" : "") + "</div>" +
					"<div class='content'>" + this.content + "</div>" +
				"</div>";

	this.display = function($container,x,y) {
		$container.append(this.html);
		var $w = $(".dragme[data-window-id='" + this.id + "']");
		$w.css("top", y);
		$w.css("left", x);
		$w.css("width", this.width);
		$w.css("height", this.height);

		$(".dragme").not($w).css("z-index", 9);
		$w.css("z-index", 10);

		if (this.close) {
			var self = this;
			$(".dragme[data-window-id='" + this.id + "'] p.x").click(function() {
				self.destroy();
			});
		}

		$w.draggable({ handle: ".handle" });

		$w.mousedown(function(e) {
				$(".dragme").not($w).css("z-index", 9);
				$w.css("z-index", 10);
		});
	};

	this.destroy = function() {	
		console.log(this.onClose);	
		if (this.onClose) {
			console.log("close callback");
			this.onClose.call();
		}
		$(".dragme[data-window-id='" + this.id + "']").remove();
	};

	this.button = function(text, action) {
		$(".dragme[data-window-id='" + this.id + "'] .handle").append("<span class='button noselect'>" + text + "</span>");
		var b = $(".dragme[data-window-id='" + this.id + "'] .handle").children(".button").last();
		b.click(action);
	};

	this.desc = function(text) {
		$(".dragme[data-window-id='" + this.id + "'] .handle").append("<span class='desc noselect'>" + text + "</span>");
	};

	this.jqObj = function() {
		return $(".dragme[data-window-id='" + this.id + "']");
	};
};
