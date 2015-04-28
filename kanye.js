/* track object constructor
 * audiotag = document.getElementById('yourAudioTag'), to be .play()'d
 * length_classes looked up in kanye.css
 * note_count = length_in_eighths/16 (2 bars)
 */
function Track (name, audiotag, lengthclass, note_count, thumbnail_url) {
	this.name = name;
	this.audiotag = audiotag;
	this.lengthclass = lengthclass;
	this.note_count = note_count;
	this.thumbnail_url = thumbnail_url;
};


function kanye_init(tracks, $ctnr) {
	var play = new Event("play");

	console.log($ctnr);
	var track;
	for (var i = 0; i < tracks.length; i++) {
		track = tracks[i];
		$ctnr.append("<form class='kanyetrack' data-name='" + track.name + "'></form>");
		var $form = $("form[data-name='" + track.name + "']");

		var note_length_in_eighths = 16 / track.note_count;

		for (var p = 0; p < track.note_count; p++) {
			$form.append("<input type='checkbox' data-track-number=" + i + " data-note-number=" + p * note_length_in_eighths + " class='kanyebox " + track.lengthclass + "' id='" + track.name + p + "'>" +
				"<label for='" + track.name + p + "'></label>");
			$form.append(track.audiotag);
			$("#" + track.name + p).on("myCustomEvent", function(event, tracks) {
				$(this).parent().children("audio").trigger("play");
			});
		}
	}
};

function kanye_play(tracks, $ctnr) {
	var n = 0;
	setInterval(function() {
		console.log($(".kanyebox[data-note-number=" + n % 16 + "]:checked"));
		$(".kanyebox[data-note-number=" + n % 16 + "]:checked").trigger("myCustomEvent", tracks);
		n++;
	}, 320); // 100 bpm
};

function kanyeWindow() {

	var width = 500;
	var height = 500;
	var kw = new Window("808s and Heartbreak Life Enhancer", "", true, width, height);
	kw.display($("body"), 100, 100);
	var $ctnr = $(".dragme[data-window-id='" + kw.id + "'] .content");

	var choralat = document.createElement("audio");
	choralat.setAttribute('src', './kanye-snips/kanye-choral-2bar.mp3');
	var choral = new Track("choral", choralat, "2bar", 1, "./kanye-snips/choral.png");

	var blocksat = document.createElement("audio");
	blocksat.setAttribute('src', './kanye-snips/kanye-blocks-half.mp3');
	var blocks = new Track("blocks", blocksat, "half", 4, "./kanye-snips/blocks.png");

	var clickat = document.createElement("audio");
	clickat.setAttribute('src', './kanye-snips/kanye-clicktrack-whole.mp3');
	var click = new Track("clicktrack", clickat, "whole", 2, "./kanye-snips/clicktrack.png");

	var growat = document.createElement("audio");
	growat.setAttribute('src', './kanye-snips/kanye-grow-2bar.mp3');
	var grow = new Track("grow", growat, "2bar", 1, "./kanye-snips/grow.png");

	var tracks = [choral, blocks, click, grow];

	kanye_init(tracks, $ctnr);
	kanye_play(tracks, $ctnr);
};