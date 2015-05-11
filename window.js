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
		$(".dragme[data-window-id='" + this.id + "'] .handle").append("<span class='button'>" + text + "</span>");
		var b = $(".dragme[data-window-id='" + this.id + "'] .handle").children(".button").last();
		b.click(action);
	};
};