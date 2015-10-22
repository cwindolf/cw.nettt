function Fourgan($ctnr) {
	console.log("fourgan!");

	$ctnr.attr("id", "fourgan")
	$ctnr.append("<ul id='sliders'></ul>");

	var $ul = $("#sliders");
	for (var i = 0; i < 9; i++) {
		$ul.append("<li><input type='range' min='0' max='1' value='0.5' class='slider' id='slider" + i + "'></li>");
		$("#slider" + i).change(update);
	}

	// set up the plot
	var width = 500; var height = 200;
	$ctnr.append("<br><canvas width='"+width+"'' height='"+height+"'' id='fourgan_canvas'>");
	var canvas = document.getElementById("fourgan_canvas");
	paper.setup(canvas);
	// add axes
	var y_axis = new paper.Path();
	y_axis.strokeColor = 'black';
	y_axis.strokeWidth = 2;
	y_axis.dashArray = [10, 4];
	var y_o = new paper.Point(width / 2, 0);
	y_axis.moveTo(y_o);
	y_axis.lineTo(new paper.Point(width / 2, height));
	var x_axis = new paper.Path();
	x_axis.strokeColor = 'black';
	x_axis.strokeWidth = 2;
	x_axis.dashArray = [10, 4];
	var x_o = new paper.Point(0, height / 2);
	x_axis.moveTo(x_o);
	x_axis.lineTo(new paper.Point(width, height / 2));

	// returns array
	function coefficients() {
		var ret = [];
		for (var j = 0; j < 9; j++) {
			ret.push(parseFloat($("#slider" + j).val()));
		}
		return ret;
	}

	// transform: [Numbers] -> (Number -> Number)
	function transform(coefts) {
		
	}

	function updatePlot(c_arr) {

		paper.view.draw();
	}

	function update() {
		var c_arr = coefficients();
		updatePlot(c_arr);
	}
	update();
}