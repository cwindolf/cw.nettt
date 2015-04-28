// x,v paper points
function PondPoint (x, v) {
	this.x = x;
	this.v = v;
}

// requires paper.js
function Pond ($ctnr) {
	this.$ctnr = $ctnr;
	this.width = Math.round($ctnr.width());
	this.height = Math.round($ctnr.height());

	this.canvas = document.createElement("canvas");
	this.canvas.setAttribute("id", "pond");
	this.canvas.setAttribute("width", this.width);
	this.canvas.setAttribute("height", this.height);
	this.$ctnr.append($(this.canvas));
	this.context = this.canvas.getContext("2d");

	paper.setup("pond");

	this.width = paper.view.width;
	this.height = paper.view.height;
	this.ripples = []; // store paths as arrays of PondPoints
	this.paths = [];

	this.renderRipple = function(ripple) {
		var path = new paper.Path();
		path.strokeColor = "black";
		for (var i = 0; i < ripple.length; i++) {
			path.add(ripple[i].x);
		}
		path.closed = true;
		path.smooth();
		this.paths.push(path);
	};

	this.render = function() {
		for (var p = 0; p < this.paths.length; p++) {
			console.log(this.paths[p]);
			this.paths[p].remove();
		}

		for (var i = 0; i < this.ripples.length; i++) {
			this.renderRipple(this.ripples[i]);
		}
	};

	this.stepRipple = function(ripple, dt) {
		var pt;
		for (var i = 0; i < ripple.length; i++) {
			pt = ripple[i];
			pt.x.add(pt.v.multiply(dt));
		}
	};

	this.step = function(dt) {
		for (var i = 0; i < this.ripples.length; i++) {
			this.stepRipple(this.ripples[i], dt);
		}
	};

	this.time = Date.now();
	var self = this;
	paper.view.on('frame', function() {
		var dt = Date.now() - self.time;
		self.step(dt);
		self.render();
		self.time = Date.now();
	});

	this.tool = new paper.Tool();
	this.tool.onMouseDown = function(event) {
		console.log("mouse down bitch");
		self.ripples.push([]);
		var pt;
		// surround event.point with path points
		for (var theta = 0; theta < Math.PI; theta += Math.PI/100) {
			pt = new paper.Point(Math.cos(theta), Math.sin(theta));
			var ppt = new PondPoint(event.point.add(pt), pt);
			self.ripples[self.ripples.length - 1].push(ppt);
		};
	}
}