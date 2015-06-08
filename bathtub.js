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
	var Q = 250;

	this.width = 650;
	this.height = 650;
	this.canvas = document.createElement("canvas");
	this.canvas.setAttribute("id", "conway");
	this.canvas.setAttribute("width", this.width);
	this.canvas.setAttribute("height", this.height);
	$ctnr.append($(this.canvas));
	this.context = this.canvas.getContext("2d");
	this.context.imageSmoothingEnabled = false;
	this.id = undefined;
	var system = new SWE(Q, this.width);
	console.log(system);
	var self = this;

	var then = Date.now();
	var now;
	var dt;

	// play sounds
	window.AudioContext = window.AudioContext||window.webkitAudioContext;
    this.acontext = new AudioContext();

	$(this.canvas).click(function(e) {
		var x = Math.floor(Math.floor((e.pageX-$(self.canvas).offset().left)) / (system.dd));
		var y = Math.floor(Math.floor((e.pageY-$(self.canvas).offset().top)) / (system.dd));
    	// pause and drip a plip at x,y
    	clearInterval(this.id);
    	system.plip(x,y);
    	then = Date.now() - 1;
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
		this.id = setInterval(function() {
			now = Date.now();
			dt = now - then;
			system.step(dt / 20000);
			then = now;
		}, 1); // fast as possible i guess. system is more stable with lower dt values.
	};

	this.cleanup = function() { clearInterval(this.id); };
	this.update();
};























