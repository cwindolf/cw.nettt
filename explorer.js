function Explorer ($ctnr) {
	this.$ctnr = $ctnr;
	this.width = Math.round(this.$ctnr.width());
	this.height = Math.round(this.$ctnr.height());
	
	// THREE.js
	this.scene = new THREE.Scene();
	this.camera = new THREE.PerspectiveCamera(75, this.width / this.height, 0.1, 1000);
	this.camera.position.set(0,0,100);
	// this.camera.lookAt(5,5,0);
	this.renderer = new THREE.WebGLRenderer({alpha:true});
	this.renderer.setSize(this.width, this.height);
	this.renderer.setClearColor( 0xffffff, 1);
	this.$ctnr.append(this.renderer.domElement);

	var light = new THREE.AmbientLight( 0x404040 ); // soft white light
	this.scene.add( light );

	//
	this.MAP_SIDE_LENGTH = 20;

	// generate terrain height data!
	this.planeG = new THREE.PlaneGeometry(60, 60, this.MAP_SIDE_LENGTH, this.MAP_SIDE_LENGTH);
	// get two functions that are MAP_SIDE_LENGTH-periodic
	this.Am = Math.floor(Math.random() * 5) + 1;
	this.An = Math.floor(Math.random() * 5) + 1;
	this.Pm = Math.floor(Math.random() * 5) + 1;
	this.Pn = Math.floor(Math.random() * 5) + 1;
	// surface as function of x,y
	this.m = function(x) {
		return this.Am * Math.cos(this.Pm * 2 * Math.PI * x / this.MAP_SIDE_LENGTH);
	};
	this.n = function(y) {
		return this.An * Math.cos(this.Pn * 2 * Math.PI * y / this.MAP_SIDE_LENGTH);
	};
	this.f = function(x, y) {
		return m(x) * n(y);
	};
	// data for every inch of the map. gonna be a huge-ass array.
	this.vertices = [];
	for (var x = 0; x < this.MAP_SIDE_LENGTH; x++) {
		var mx = this.m(x);
		for (var y = 0; y < this.MAP_SIDE_LENGTH; y++) {
			this.planeG.vertices[x * this.MAP_SIDE_LENGTH + y].z = Math.floor(mx * this.n(y));
		}
	}

	this.planeM = new THREE.MeshLambertMaterial({
		color: 0x8899aa,
		// wireframe: true
	});
	this.planeM.shading = THREE.FlatShading;
	this.plane = new THREE.Mesh(this.planeG, this.planeM);
	this.scene.add(this.plane);

	var self = this;
	function render() {
		requestAnimationFrame( render );
		self.plane.rotation.y += .01;
		// do stuff
		self.renderer.render( self.scene, self.camera );
	}
	render();
}