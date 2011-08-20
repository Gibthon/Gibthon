
function px2em(input) 
{
    var emSize = parseFloat($('body').css("font-size"));
    return (input / emSize);
}

function em2px(input) 
{
    var emSize = parseFloat($('body').css("font-size"));
    return (input * emSize);
}


var meta_initial_html = '' +
'<div class="ui-widget-content ui-corner-top top-box" >' +
'	<h3 id="desc"></h3>' +
'	<h5 id="origin"></h5>' +
'</div>' +
'<div class="ui-widget-content ui-corner-bottom bottom-box" >' +
'	<span id="more_details">' +
'		<span id="more">Show</span> Details ' +
'		<span id="icon" style="display:inline-block;padding-top:.2em;height:0.8em;" class="ui-icon ui-icon-triangle-1-s" />' +
'	</span>' +
'	<div id="drop_down">' +
'		<h4>Annotations</h4>' +
'			<div class="info" id="annot_div">' +
'				<table id="annot_table">' +
'				</table>' +
'			</div>' +
'		<h4>References</h4>' +
'			<div class="info" id="ref_div">' +
'				' +
'			</div>' +
'	</div>' +
'</div>';

(function( $, undefined ) {

$.widget("ui.fragmentMeta", {
	options: {
		id: 0,
	},
	_create: function() {
		//Init the element, and fetch the initial data
		var self = this;
		this.$el = $(this.element[0]).html(meta_initial_html);
		this.$dd = this.$el.find('#drop_down');
		this.$more_btn = this.$el.find('#more_details')
			.click(function() {self.show();});
		this.$more = this.$el.find('#more');
		this.$desc = this.$el.find('#desc');
		this.$origin = this.$el.find('#origin');
		this.$icon = this.$el.find('#icon');
		
		this.$annotation = this.$el.find('#annot_table');
		this.$ref = this.$el.find('#ref_div');
		
		//fetch name and origin
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'meta',}, function(data) {
			if(data[0] == 0)
			{
				self.$desc.text(data[1].desc);
				self.$origin.text(data[1].length + "bp, " + data[1].origin);
			}
			else
			{
				console.log(data[1] + ", while getting meta.");
			}
		});
		
		//fetch annotations
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'annotations',}, function(data) {
			if(data[0] == 0)
			{
				var odd = false;
				for (key in data[1])
				{
					self.$annotation.append(self._make_annotation(key, data[1][key], odd));
					odd = !odd;
				}
			}
			else
			{
				console.log(data[1] + ", while getting annotations.");
			}
		});
		
		//fetch references
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'refs',}, function(data) {
			if(data[0] == 0)
			{
				for (key in data[1])
				{
					self.$ref.append(self._make_reference(data[1][key]));
				}
				$('.ref-search-btn').button({
														icons: {
															primary: 'ui-icon-search'
																},
														label: 'Search',
													});
			}
			else
			{
				console.log(data[1] + ", while getting references.");
			}
		});
		
	},
	show: function() {
		var self = this;
		this.$more.text('Hide');
		this.$more_btn
			.unbind()
			.click(function() {self.hide();});
		this.$icon
			.removeClass('ui-icon-triangle-1-s')
			.addClass('ui-icon-triangle-1-n');
		this.$dd.slideDown(250);
	},
	hide: function() {
		var self = this;
		this.$more.text('Show');
		this.$more_btn
			.unbind()
			.click(function() {self.show();});
		this.$icon
			.removeClass('ui-icon-triangle-1-n')
			.addClass('ui-icon-triangle-1-s');
		this.$dd.slideUp(250);
	},
	_make_annotation: function(key, value_list, alt) {
		var cls = 'tr';
		if(alt == true)
			cls = 'tr-alt';
		var ret = "<tr class='" + cls + "'><td>" + key + ": </td><td>";
		for (i in value_list)
		{
			ret = ret + value_list[i];
			if(i != value_list.length -1)
				ret = ret + ", ";
		}
		return ret + "</td></tr>";
	},
	_make_reference: function(ref)
	{
		var ret = '' +
		'<div class="ref">' +
		'	<div class="ref-title-wrap">' +
		'		<div class="ref-title-div">' +
		'			<h4 class="ref-title">' +
		'				' + ref.title +
		'			</h4>' +
		'		</div>' +
		'		<div class="ref-search-div">' +
		'			<button class="ref-search-btn" onclick="window.open(\'http://www.google.com?q=' + ref.title + '\',\'_blank\');"></button>' +
		'		</div>' +
		'		<div style="clear:both;"></div>' +
		'	</div>' +
		'	<div class="ref-detail">' +
		'		<div class="ref-authors-div">' +
		'			<h6 class="ref-authors">' + ref.authors + '</h6>'
		'		</div>' +
		'		<div class="ref-journal-div">' +
		'			<h6 class="ref-journal">' + ref.journal + '</h6>' +
		'		</div>' +
		'		<div class="ref-links">';
		
		if(ref.medline_id != "")
		{
			ret = ret + '<h6 class="ref-link">Meline ID: ' + ref.medline_id + '</h6>';
		}
		if(ref.pubmed_id != "")
		{
			ret = ret + '<h6 class="ref-link">PubMed ID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + 
					ref.pubmed_id + '" target="_blank">' + ref.pubmed_id + '</a></h6>';
		}
		ret = ret + 
		'		</div>' +
		'	</div>' +
		'</div>';
		return ret;
	},
});

})( jQuery );	

var initial_sequence_html = '' +
'<div id="loader" class="ui-state-highlight content middle-box">' +
'	<span id="load_span"><b>Loading:</b> <span id="progress">0</span>/<span id="length">0</span> bp </span>' +
'	<span id="progress_span">' +
'		<div id="progressbar" ></div>' +			
'	</span>' +
'</div>' +
'<div class="ui-widget-content ui-corner-bottom bottom-box content">' +	
'	<div id="seq_wrap unselectable" unselectable="on">' +
'		<div id="seq_inner" class="unselectable" unselectable="on">' +
'		</div>' +
'	</div>' +
'</div>';

//inital guess at char width
var char_width = 0.64;
var click_thresh = 300; //ms

var select_start = '<span class="selected unselectable" unselectable="on">';
var select_end = '</span>';

(function( $, undefined ) {

$.widget("ui.fragmentSequence", {
	options: {
		id: 0,
	},
	_create: function() {
		var self = this;
		this.$el = $(this.element[0]).html(initial_sequence_html);
		this.$len = this.$el.find('#length');
		this.$prog = this.$el.find('#progress');
		this.$bar = this.$el.find('#progressbar').progressbar({value: 0,});
		this.$barval = this.$el.find('.ui-progressbar-value');
		this.$loader = this.$el.find('#loader');
		this.$seq = this.$el.find('#seq_inner');
		this.len = 0;
		this.seq = "";
		this.pos = 0;
		this.rowlength = 10 * parseInt(px2em(this.$el.width()) / (11 * char_width));
		this._get_seq_meta();
		
		this._get_char_width();
		this._mouse_down = false;
		
		$(window).mouseup(function(event) {self._window_mouse_up();});
				
		this.rows = new Array();
	},
	_get_seq_meta: function(){
		var self = this;
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'seq_meta',}, function(data) {
			if(data[0] != 0)
			{
				console.log(data[1] + " while getting sequence metadata");
				return;
			}
			self.len = data[1].len;
			self.$len.text(data[1].len);
			if(self.len > 2000) //if we should load progressively 
			{
				self.$loader.slideDown(100);
			}
			
			self.features = data[1].feats;
			self.alphabet = data[1].alpha;
			self.alphabet[' '] = ' ';
			
			self._get_seq(0);
		});
	},
	_get_seq: function(offset){
		var self = this;
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'seq', 'offset': offset}, function(data) {
			if(data[0] != 0)
			{
				console.log(data[1] + " while getting seq");
				return;
			}
			self.seq = self.seq + data[1];
			var flush = false;//should we flush all the sequence?
			if(self.seq.length < self.len)//is there more data to come?
			{
				//get some more data
				self._get_seq(self.seq.length);
			}
			else
			{
				flush = true;
				self.$loader.slideUp(500);
			}
			
			//while we have enough data to make a complete row
			while((self.seq.length - self.pos) > self.rowlength)
			{
				self.pos = self.pos + self._make_row(self.pos);
			} 
			if(flush)
			{
				self.pos = self.pos + self._make_row(self.pos);
				self._label_features();
				self._get_char_width();
				//for(i in self.rows)
				//{
				//	console.log('this.rows[' + i + '].find(".seq-fwd")text() = ' + self.rows[i].find('.seq-fwd').text());
				//}
			}
			
			self.$prog.text(self.seq.length);
			self.$bar.progressbar('value', parseInt((100 * self.seq.length) / self.len));
		});
	},
	_complement: function(string){
		var ret = "";
		for(i in string)
		{
			ret = ret + this.alphabet[string[i]];
		}
		return ret;
	},
	_make_row: function(start){
		var self = this;
		var s = "";
		var end = start + this.rowlength;
		if( start == 0)
		{
			s = " " + this.seq.substr(start, this.rowlength - 1);
			end = end - 1;
		}
		else
			s = this.seq.substr(start - 1, this.rowlength);
		
		var seq = "";
		var label = "";
		for(var i = 0; i < this.rowlength; i = i + 5)
		{
			seq = seq + s.substr(i, 5) + " ";
			label = label + '<div class="seq-label">' + (start + i) + '</div>';
		}
		var cseq = this._complement(seq);
		
		var fwd_feats = "";
		var rev_feats = "";
		for(f in this.features)
		{
			var feat_html = make_feat_html(start, start + this.rowlength, this.features[f], f);

			if(this.features[f].strand > 0) fwd_feats = fwd_feats + feat_html;
			else rev_feats = rev_feats + feat_html;
			
		}
		
		
		var $row = $(document.createElement('div'))
			.attr('unselectable', 'on')
			.addClass('row unselectable')
			.html(	'<div class="ladder unselectable" unselectable="on"></div>' +
					'<div id="' + start + '-fwd-feat" class="feat-ref feat unselectable" unselectable="on">' + fwd_feats + '</div>' +
					'<div style="position:relative;">'  +
					'	<div id="' + start + '-fwd-hid" class="seq-fwd-hid seq-hid unselectable" unselectable="on">' + seq + '</div>' +
					'	<div id="' + start + '-fwd" class="seq-fwd seq unselectable" unselectable="on">' + seq + '</div>' +
					'	<div id="' + start + '-label" class="label-div unselectable" unselectable="on">' + label + '</div>' +
					'	<div id="' + start + '-rev" class="seq-rev seq unselectable" unselectable="on">' + cseq + '</div>' +
					'</div>' +
					'<div id="' + start + '-rev-feat" class="feat-rev feat unselectable" unselectable="on">' + rev_feats + '</div>' +
					'<div class="ladder unselectable" unselectable="on"></div>'
			);
		
		$row.find('.seq-fwd')
			.mousedown	(function(event) {self._on_mouse_down(event);})
			.mouseup 	(function(event) {self._on_mouse_up	 (event);})
			.mousemove	(function(event) {self._on_mouse_move(event);});
		
		this.$seq.append($row);
		this.rows.push($row);
	//console.log('pushed to rows, length: ' + this.rows.length + ', $row.find(".seq-fwd")text() = ' + $row.find('.seq-fwd').text());
		return this.rowlength;
	},
	_label_features: function()
	{
		for(f in this.features)
		{
			$('.feat-id-' + f).tipTip({
				maxWidth: '800px',
				content: make_feat_tooltip(this.features[f]),
			});
		}
	},
//SELECTION FUNCTIONS
	_select: function(start, end)
	{
		if(start > end)
		{
			var t = end; end = start; start = t;
		}
		if(this._selected && (this._select_start == start) && (this._select_end == end))
			return;
		
		this._select_start = start;
		this._select_end = end;
				
		this._selected = true;		
		//work out the start and end rows
		var sr = parseInt(Math.floor(start / this.rowlength));
		var er = parseInt(Math.floor(end / this.rowlength));
		
		var $sr = this.rows[sr].contents().children('.seq-fwd');
		var $srh = this.rows[sr].contents().children('.seq-fwd-hid');
		var $er = this.rows[er].contents().children('.seq-fwd');
		var $erh = this.rows[er].contents().children('.seq-fwd-hid');
		
		sr_ = toPadded(start - this.rowlength * sr);
		er_ = toPadded(end   - this.rowlength * er);
		
		//sr == er
		if(sr == er)
		{
 			var html = $srh.text().substr(0, sr_) + 
						select_start + 
							$srh.text().substr(sr_, er_ - sr_) + 
						select_end +
						$srh.text().substr(er_);
			$sr.html(html);
			return;
		}
		
		//change start row
		$sr.html( $srh.text().substr(0, sr_) + select_start + $srh.text().substr(sr_) + select_end);
		
		//change imtemediate rows
		for(var row = sr + 1; row < er; row = row + 1)
		{
			this._select_row(row);
		}
		
		//change end row
		$er.html( select_start + $erh.text().substr(0, er_) + select_end + $erh.text().substr(er_));
	},
	_select_row: function(row)
	{
		var $row = this.rows[row].contents().children('.seq-fwd');
		var $rowh= this.rows[row].contents().children('.seq-fwd-hid');
		$row.html(select_start + $rowh.text() + select_end);
	},
	_clear_selection: function()
	{
		this._selected = false;
		$('.selected').parents('.seq-fwd').each( function()
		{
			$(this).html( $(this).prev().text());
		});
	},
//MOUSE FUNCTIONS
	_on_mouse_down: function(event)
	{
		this._mouse_down_stamp = event.timeStamp;
		this._mouse_down_pos = this._get_mouse_pos(event);
		this._mouse_down = true;
	},
	_on_mouse_move: function(event)
	{
		if(this._mouse_down)
			this._select(this._mouse_down_pos, this._get_mouse_pos(event));
	},
	_on_mouse_up: function(event)
	{
		this._mouse_down = false;
		if((event.timeStamp - this._mouse_down_stamp) < click_thresh)
		{
			//treat the event as a click
			if(this._selected)
			{
				this._clear_selection();
				return;
			}
			else
			{
				//select the clicked block
				var pos = this._get_mouse_pos(event);
				start = 5 * Math.floor(pos/ 5);
				this._select(start, start + 5);
				return;
			}
		}
		var pos = this._get_mouse_pos(event);
		if(this._selected) this._clear_selection();
		this._select(this._mouse_down_pos, pos);
	},
	_window_mouse_up: function()
	{
		this._mouse_down = false;
	},
	_get_mouse_pos: function(event) //given a mouse event, return the offset of the clicked character.
	{
		var largenum = parseInt(event.currentTarget.id.split('-')[0]);
		
		var smallnum = Math.floor((event.pageX - $(event.target).offset().left) / this.char_width);
		smallnum = toUnPadded(smallnum);
		
		//if($(event.target).hasClass('selected')) console.log('selected position:' + largenum + ' + ' + smallnum + " event.currentTarget.id=" + event.currentTarget.id);
		return largenum + smallnum;
	},
	_get_char_width: function()
	{
		$r = $('#0-fwd');
		
		if($r.length == 0) this.char_width = char_width;
		else
		{
			this.char_width = $r.width() / $r.text().length;
		}
		return this.char_width;
	}
});

})( jQuery );	


/*
*
**** Helper functions
*
*/
var make_feat_html = function(r_s, r_e, feature, f_id)
{
	var f_s = feature.start + 1; var f_e = feature.end + 1; //features stored as 0-offset
	if(f_e < f_s) //make certain the end and the start are the right way around
	{
		var t = f_e; f_e = f_s; f_s = t;
	}
	//check if the feature is not present in the row
	if((f_e < r_s) || (f_s > r_e))
	{
		return "";
	}
	
	var left = f_s - r_s;
	if(left < 0) left = 0;
	var right = r_e - f_e;
	if(right < 0) right = 0;
	
	var l = r_e - r_s;
	
	var r = "";
	for(var i = 0; i < l; i = i + 5)
	{
		for(var j = 0; j < 5; j = j + 1)
		{
			if((i+j) == left)
			{
				r = r + '<span class="feat-hl feat-type-' + feature.type + ' feat-id-' + f_id + '">';
			}
			if((l - (i+j)) == right)
			{
				r = r + '</span>';
			}
			r = r + ' ';
		}
		
		if((i+j) < l)
			r = r + ' ';
	}
	
	return '<div class="feat-div">' + r + '</div>';	
}
var make_feat_tooltip = function(feature)
{
	var strand = '';
	if(feature.strand == 1)
		strand = 'fwd';
	else if(feature.strand == -1)
		strand = 'rev';
		
	var ret = '' +
	'<table style="word-wrap:break-word;">' +
	'	<tr><td style="min-width:7em;">Type :</td><td>' + feature.type + '</td></tr>' +
	'	<tr><td>Location :</td><td>' + (feature.start) + '-' + (feature.end) + '</td></tr>' +
	'	<tr><td>Length</td><td>' + (feature.end - feature.start) + '</td></tr>' +
	'	<tr><td>Strand :</td><td>' + strand + '</td></tr>';
	
	for(i in feature.qualifiers)
	{
		ret = ret + '<tr><td>' + feature.qualifiers[i].name + '</td><td>' + feature.qualifiers[i].data + '</td></tr>';
	}
	ret = ret + '</table>';
	return ret;
}
var toPadded = function(len)
{
	return 6 * Math.floor(len / 5) + (len % 5);
}
var toUnPadded = function(len)
{
	return 5 * Math.floor(len / 6) + (len % 6);
}
