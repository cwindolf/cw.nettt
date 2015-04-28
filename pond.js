// x,v paper points
function PondPoint (x, v) {
	this.x = x;
	this.v = v;
}

// requires paper.js
function Pond ($ctnr) {
	var self = this; // <3 JS!
	this.INSET = 10;
	this.E = .8;
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

	this.ripples = []; // store paths as arrays of PondPoints
	this.pathsGroup = new paper.Group();

	var spoint = new paper.Point(this.INSET,this.INSET);
	var ssize = new paper.Size(self.width - this.INSET * 2, self.height - this.INSET * 2);;
	this.shore = new paper.Path.Rectangle(spoint, ssize);
	this.shore.strokeColor = 'black';

	this.render = function() {
		this.pathsGroup.removeChildren();

		for (var i = 0; i < this.ripples.length; i++) {
			var path = new paper.Path();
			path.strokeColor = "black";
			for (var p = 0; p < this.ripples[i].length; p++) {
				path.add(this.ripples[i][p].x);
			}
			path.closed = true;
			path.smooth();
			this.pathsGroup.addChild(path);
		}
	};

	this.step = function(dt) {
		for (var i = 0; i < this.ripples.length; i++) {
			var pt;
			for (var p = 0; p < this.ripples[i].length; p++) {
				pt = this.ripples[i][p];
				if (!this.shore.contains(pt.x)) {
					// flip y-velocity
					if (pt.x.y <= this.INSET || pt.x.y >= this.height - this.INSET * 2) {
						pt.v.y = -this.E * pt.v.y;
						pt.x.y += 0.1 * pt.v.y;
					}
					// flip x-velocity
					if (pt.x.x <= this.INSET || pt.x.x >= this.width - this.INSET * 2) {
						pt.v.x = -this.E * pt.v.x;
						pt.x.x += 0.1 * pt.v.x;
					}
				}
				pt.x = pt.x.add(pt.v.multiply(dt / 1000)); //ms
			}
		}
	};

	this.time = Date.now();
	paper.view.on('frame', function() {
		if (self.ripples.length > 0) {
			var dt = Date.now() - self.time;
			self.step(dt);
			self.render();
			self.time = Date.now();
		}
	});

	this.tool = new paper.Tool();
	this.tool.onMouseDown = function(event) {
		self.time = Date.now();
		self.ripples.push([]);
		var pt;
		// surround event.point with path points
		for (var theta = 0; theta < 2* Math.PI; theta += Math.PI/128) {
			pt = new paper.Point(Math.cos(theta), Math.sin(theta));
			var ppt = new PondPoint(event.point.add(pt), pt.multiply(30));
			self.ripples[self.ripples.length - 1].push(ppt);
		};
	}
}