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