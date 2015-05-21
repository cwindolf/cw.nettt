function Explorer ($ctnr) {
	this.$ctnr = $ctnr;
	this.width = Math.round(this.$ctnr.width());
	this.height = Math.round(this.$ctnr.height());
	
	// THREE.js
	this.scene = new THREE.Scene();
	this.camera = new THREE.PerspectiveCamera(75, this.width / this.height, 0.1, 1000);
	this.camera.position.set(15,-950,150);
	this.camera.rotation.x += 1.5; 

	// this.camera.lookAt(5,5,0);
	this.renderer = new THREE.WebGLRenderer({
		alpha:true,
		precision: "highp",
	});
	this.renderer.setSize(this.width, this.height);
	this.renderer.setClearColor( 0xbbbbbb, 1);
	this.$ctnr.append(this.renderer.domElement);

	var alight = new THREE.AmbientLight( 0x404040 ); // soft white light
	alight.castShadow = true;
	alight.intensity = 0.5;
	this.scene.add( alight );

	this.plight = new THREE.PointLight( 0xffddee, 0.9, 300 );
	this.plight.castShadow = true;
	this.plight.position.set(15,-950, 160);
	this.scene.add(this.plight);

	//
	this.MAP_SIDE_LENGTH = 500;

	// generate terrain height data!
	this.planeG = new THREE.PlaneGeometry(1000, 1000, this.MAP_SIDE_LENGTH, this.MAP_SIDE_LENGTH);
	// get two functions that are MAP_SIDE_LENGTH-periodic
	this.Am = 15 + 1;
	this.An = 5 + 1;
	this.Pm = Math.floor(Math.random() * 20) + 1;
	this.Pn = Math.floor(Math.random() * 30) + 1;
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
	var vertices = new Float32Array( this.MAP_SIDE_LENGTH * this.MAP_SIDE_LENGTH );
	for (var x = 0; x < this.MAP_SIDE_LENGTH; x++) {
		var mx = this.m(x + 1);
		for (var y = 0; y < this.MAP_SIDE_LENGTH; y++) {
			this.planeG.vertices[x * this.MAP_SIDE_LENGTH + y].z = Math.floor(mx * this.n(y + 1));
		}
	}
	this.planeG.mergeVertices();
	this.planeG.computeFaceNormals();
	this.planeG.computeVertexNormals();
	this.planeM = new THREE.MeshLambertMaterial({
		color: 0x0000ff,
		// wireframe: true
	});
	this.plane = new THREE.Mesh(this.planeG, this.planeM);
	this.plane.castShadow = true;
	this.plane.receiveShadow = true;
	this.scene.add(this.plane);

	var self = this;
	function render() {
		requestAnimationFrame( render );
		// self.plane.rotation.y += .01;
		self.camera.position.y += 1;
		self.plight.position.y += 1;
		// do stuff
		self.renderer.render( self.scene, self.camera );
	}
	render();
}