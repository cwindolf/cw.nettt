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