
var httppat = /^http:\/\//i

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

(function( $, undefined ) {

$.widget("ui.fragmentMeta", {
	options: {
		frag: undefined,
	},
	_create: function() {
		//Init the element, and fetch the initial data
		var self = this;
		this.$el = $(this.element[0]);
		this.$dd = this.$el.find('#drop_down');
		this.$more_btn = this.$el.find('#more_details')
			.click(function() {self.show();});
		this.$more = this.$el.find('#more');
		this.$name = this.$el.find('#name').children('.magic-text');
		this.$desc = this.$el.find('#desc').children('.magic-text');
		this.$origin = this.$el.find('#origin');
		this.$icon = this.$el.find('#icon');
		
		this.$annotation = this.$el.find('#annot_table');
		this.$annot_div = this.$el.find('#annot_div');
		this.$form = this.$el.find('form');
		this.$ref = this.$el.find('#ref_div');
		
		self.$form.magicForm({
			edit: function() {
				self.was_visible = self.visible;
				self.show();
				self.$annot_div.formExtender('enable');
			},
			cancel: function() {
				if(!self.was_visible) self.hide();
				self.$annot_div.formExtender('disable');
			},
			submit_on_save: false,
			save: function() {
				if(!self.was_visible) self.hide();
				self.$annot_div.formExtender('disable');
				
				var meta = self.getMeta();
				self.fragment.setMeta(meta, function(){}, function(err) {console.error(
					'Error setting metadata "'+err+'"');});
			},
		});
		self.$annot_div.formExtender({
			unique: true,
			disabled: true,
			addInitial: false,
			beforeAdd: function(e, $el)
			{
				$el.find('.magic-item').each( function() {self.$form.magicForm('SetState', $(this));});
			},
			add: function(e, $el) {
				self._remark();
			},
			remove: function(e, $el) {
				self._remark();
			},
		});
		
		this.fragment = this.options.frag;
		this.fragment.getMeta( function(meta)
							 {
								 self._display_metadata(meta);
							 });

		this.visible = false;
		this.was_visible = false;			
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
		this.$dd.slideDown(250, function(){self.visible=true;});
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
		this.$dd.slideUp(250, function(){self.visible=false;});
	},
	_remark: function()
	{
		var self = this;
		var odd = false;
		this.$annotation.find('tr').each( function() {
			$(this).removeClass('tr tr-alt');
			if(odd) $(this).addClass('tr-alt')
			else $(this).addClass('tr');
			odd = !odd;
		});
	},
	_display_metadata: function(meta) {
		var self = this;
		
		//clear any previous metadata
		self.$annotation.html('');
		self.$ref.html('')
		
		//Display name & description
		self.$name.text(meta.name);
		self.$desc.text(meta.desc);
		self.$origin.text(meta.length + "bp, " + meta.origin);
				
		//Display Annotations
		var odd = false;
		for (key in meta.annots)
		{
			self.$annot_div.formExtender('Add', [key, meta.annots[key]]);
			odd = !odd;
		}
		
		//Display References
		for (r in meta.refs)
		{
			self.$ref.append(self._make_reference(meta.refs[r]));
		}
		$('.ref-search-btn').button({
			icons: {primary: 'ui-icon-search'},
			label: 'Search',
		});
	},
	_make_annotation: function(key, v, alt) {
		var cls = 'tr';
		if(alt == true)
			cls = 'tr-alt';
		var value = '';
		if(v.constructor == Array) //if v is an array
		{
			for (i in value_list)
			{
				value = value + value_list[i];
				if(i != value_list.length -1)
					value = value + ", ";
			}
		}
		else
			value = v;
		
		if(httppat.test(value))
		{
			value = "<a href='" + value + "'>" + value + "</a>";
		}
		
		var ret = "" +
		"<tr class='" + cls + " extender-item'>" +
		"	<td id='annot_key' class='magic-item annot-key' >" +
		"		<span class='magic-text'>" + key + "</span>" +
		"		<input name='annot_key' class='magic-input' style='display:none;' type='text' value='Key' size=10 ></input>" +
		"	</td>" +
		"	<td id='annot_value' class='magic-item'>" +
		"		<span class='magic-text'>" + value + "</span>" +
		"	</td>" +
		"	<td><button class='extender-remove'>Remove</button></td>"
		"</tr>";
		return ret;
	},
	getAnnots: function()
	{
		var self = this;
		var ret = [];
		
		$.each(self.$annotation.find('tr'), function(i, val){
			var v = $(val);
			ret.push([
				v.find('#annot_key span').text(), 
				v.find('#annot_value span').text()
			]);
		});
		return ret;
	},
	_make_reference: function(ref)
	{
		var ret = '' +
		'<div class="ref">' +
		'	<div class="ref-title-wrap">' +
		'		<div class="ref-title-div">' +
		'			<h4 class="ref-title">' + ref.title + '</h4>' +
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
			ret = ret + '<h6 class="ref-mline ref-link">Meline ID: ' + ref.medline_id + '</h6>';
		}
		if(ref.pubmed_id != "")
		{
			ret = ret + '<h6 class="ref-pmed ref-link">PubMed ID: <a href="http://www.ncbi.nlm.nih.gov/pubmed/' + 
					ref.pubmed_id + '" target="_blank">' + ref.pubmed_id + '</a></h6>';
		}
		ret = ret + 
		'		</div>' +
		'	</div>' +
		'</div>';
		return ret;
	},
	getRefs: function()
	{
		var self = this;
		var ret = [];
		
		$.each(self.$ref.find('div.ref'), function(i, val){
			var v = $(val);
			ret.push( 
				new Reference(
					v.find('.ref-title').text(), 
					v.find('.ref-authors').text(),
					v.find('.ref-journal').text(),
					v.find('.ref-mline').text(),
					v.find('.ref-pmed').text()
				)
			);
		});
		return ret;
	},
	getMeta: function() {
		var self = this;
		//get metadata
		var meta = {
			'id': self.options.id, 
			'name': self.$name.text(), 
			'desc': self.$desc.text(), 
			'refs': self.getRefs(), 
			'annots': self.getAnnots(),
		};
		return meta;
	},
});

})( jQuery );	

//inital guess at char width
var char_width = 0.55;
var click_thresh = 300; //ms

var select_start = '<span class="selected unselectable" unselectable="on">';
var select_end = '</span>';

(function( $, undefined ) {

$.widget("ui.fragmentSequence", {
	options: {
		frag: undefined,
	},
	_create: function() {
		var self = this;
		this.$el = $(this.element[0]);
		this.fragment = this.options.frag;
		if(this.fragment==undefined)
			console.error('ui.fragmentSequence: Handed undefined Fragment object');
		
		this.$el.find('#copy_btn').button({
			label: "Copy Selection",
			icons: {primary: 'ui-icon-copy'}
		}).click(function() {
			show_copy_dialog();
		});
		this.$el.find('#select_btn').button({
			label: "Select All",
			icons: {primary: 'ui-icon-triangle-2-e-w'}
		}).click( function() {
			self._select(1, self.len + 1);
		});
        this.$el.find('#crop_btn').button({
            label: 'Crop to Selection',
            icons: {primary: 'ui-icon-scissors',},
        }).click( function() {
            self.cropform.cropform('setRange', self._select_start, 
                                  self._select_end);
            $('#crop_dlg').dialog('open');
        });
		
		this.$el.find('#view').buttonset();
		this.$el.find('#ds').click(function(){
			self.$el.find('.seq-fwd').removeClass('seq-hid');
			self.$el.find('.seq-rev').removeClass('seq-hid');
		});
		this.$el.find('#ss').click(function(){
			self.$el.find('.seq-fwd').removeClass('seq-hid');
			self.$el.find('.seq-rev').addClass('seq-hid');
		});
		this.$el.find('#ns').click(function(){
			self.$el.find('.seq-fwd').addClass('seq-hid');
			self.$el.find('.seq-rev').addClass('seq-hid');
		});

		//this._toolbar_height = $('#seq_toolbar').offset().top;
		this.toolbar_fixed = false;
				
		this.$len = this.$el.find('#length');
		this.$prog = this.$el.find('#progress');
		this.$bar = this.$el.find('#progressbar').progressbar({value: 0,});
		this.$barval = this.$el.find('.ui-progressbar-value');
		this.$loader = this.$el.find('#loader');
		this.$seq = this.$el.find('#seq_inner');
		this.len = 0;
		this.seq = "";
		this.pos = 0;
		this.rowlength = parseInt(px2em(this.$el.width()) / char_width);
		this.rowlength = 5 * Math.floor(this.rowlength / 6);


		this._get_seq_meta();
		
		this._get_char_width();
		this._mouse_down = false;
		
		$(window)
			.mouseup(function(event) {self._window_mouse_up();})
			.scroll(function(event) {self._on_scroll();});
				
		this.seq_disp = new Array();
		this.seq_text = new Array();
		this._selected = false;
		this._update_toolbar();

        $('#crop_dlg').dialog({
            title: 'Crop',
            autoOpen: false,
            modal: true,
            draggable: false,
            resizable: false,
            width: 450,
        });
	
        this.cropform = $('#cropform').cropform({
            fragment: this.fragment,
            cancel: function()
            {
                $('#crop_dlg').dialog('close');
            },
        });
	
	},
	_get_seq_meta: function(){
		var self = this;
		this.fragment.getFeats(function(d) 
		{
			self.len = d.length;
			self.$len.text(d.length);
			if(self.len > 1024) //if we should load progressively 
			{
				self.$loader.show();
			}
			
			self.features = d.feats;
			
			self._get_seq();
		});
	},
	_get_seq: function()
	{
		var self = this;
		if(this.fragment == undefined)
		{
			console.error('Cannot get sequence for undefined Fragment');
			return;
		}

		this.fragment.getSequence(function(seq) //update function
		{
			self.seq = ' ' + seq;
		
			//while we have enough data to make a complete row
			while((self.seq.length - self.pos) > self.rowlength)
            {
				self.pos = self.pos + self._make_row(self.pos);
            }
		
			self.$prog.text(self.seq.length);
			self.$bar.progressbar('value', parseInt((100 * self.seq.length) / self.len));
		}, function(seq) //Complete function
		{
			self.seq = ' ' + seq;
			//we're done, so remove the progress bar
			self.$loader.slideUp(500);
			
			//Use up all the remaining data
			while(self.seq.length > self.pos)
            {
				self.pos = self.pos + self._make_row(self.pos);
            }

			self._label_features();
			self._get_char_width();
            self.$seq.width(toPadded(self.rowlength) * self.char_width);
		});

	},
	_make_row: function(start){
		var self = this;
        var len = this.rowlength;
		var end = start + this.rowlength;

		var s = this.seq.substr(start, len);
		
		var seq = "";
		var label = "";
		for(var i = 0; i < s.length; i = i + 5)
		{
			seq = seq + s.substr(i, 5) + " ";
			var left= '';
			for(var j = 0; j < toPadded(i); j = j+1)
				left = left + ' ';
			label = label + '<div class="seq-label unselectable" unselectable="on" style="left:0;">' 
				+ left + '<span class="seq-label-text">' + (start + i) + '</span></div>';
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
					'<div class="sequenceStrands" id="' + start + '-seq-strands" style="position:relative;">'  +
					'	<div id="' + start + '-fwd" class="seq-fwd seq unselectable" unselectable="on">' + seq + '</div>' +
					'	<div id="' + start + '-label" class="label-div unselectable" unselectable="on">' + label + '</div>' +
					'	<div id="' + start + '-rev" class="seq-rev seq unselectable" unselectable="on">' + cseq + '</div>' +
					'</div>' +
					'<div id="' + start + '-rev-feat" class="feat-rev feat unselectable" unselectable="on">' + rev_feats + '</div>' +
					'<div class="ladder unselectable" unselectable="on"></div>'
			);
		
		$row.find('.sequenceStrands')
			.mousedown	(function(event) {self._on_mouse_down(event);})
			.mouseup 	(function(event) {self._on_mouse_up	 (event);})
			.mousemove	(function(event) {self._on_mouse_move(event);});
		
		this.$seq.append($row);
		
		this.seq_disp.push($row.find('.seq-fwd'));
		this.seq_text.push(seq);
	
		return len;
	},
	_label_features: function()
	{
		var self = this;
		for(var f in this.features)
		{
			var $feat = $('.feat-id-' + f);
			$feat.tipTip({
				maxWidth: '800px',
				content: make_feat_tooltip(this.features[f]),
			})
			.hover(function(event) { //in
				var fid = $(event.target).attr('id').split('-')[0];
				$('.feat-id-' + fid).addClass('feat-hover');
			}, 
			function(event) {//out
				var fid = $(event.target).attr('id').split('-')[0];
				$('.feat-id-' + fid).removeClass('feat-hover');
			})
			.click(function(event) {
				var fid = $(event.target).attr('id').split('-')[0];
				var feat = self.features[fid];
				self._select(feat.start + 1, feat.end + 1);
			});
		}
	},
//SELECTION FUNCTIONS
	_clear_selection: function()
	{
		if(this._selected)
		{
			$rm = this.$el.find('.selected');
			if($rm.length > 0)
			{
			var self = this;
			$rm.each( function(i, item) {
				row =  parseInt($(this).parent().attr('id').split('-')[0]);
				row = Math.floor( row / self.rowlength);
				
				self.seq_disp[row].html(self.seq_text[row]);
			});
		}
		}
		this._selected = false;
		this._update_toolbar();
	},
	_select: function(start, end)
	{
		if(!this._selected)
		{
			this._selected = true;
			this._update_toolbar();
		}
		if(start > end) {var t = end; end = start; start = t;}
		if((start == this._select_start) && (end == this._select_end))
			return;
		
		
		this._select_start = start;
		this._select_end = end;
		
		var sr = Math.floor(start / this.rowlength);
		var er = Math.floor(end / this.rowlength);
		
		//tag old selection for removal
		this.$el.find('.selected').addClass('oldSelection');
		
		if(sr == er) //special case
		{
			this._select_row(sr, start, end);
		}
		else
		{
			this._select_row(sr, start);
			for(var i = sr + 1; i < er; i = i + 1)
				this._select_row(i);
			this._select_row(er, 0, end);
		}
		
		$rm = this.$el.find('.oldSelection');
		if($rm.length > 0)
		{
			var self = this;
			$rm.each( function(i, item) {
				row =  parseInt($(this).parent().attr('id').split('-')[0]);
				row = Math.floor( row / self.rowlength);
				
				self.seq_disp[row].html(self.seq_text[row]);
			});
		}
	},
	_select_row: function(row, start, end)
	{
		if(start == undefined)
			start = 0;
		else
			start = toPadded(start % this.rowlength);
		if(end == undefined)
			end = toPadded(this.rowlength);
		else
			end = toPadded(end % this.rowlength);
		
		var t = this.seq_text[row];
		this.seq_disp[row].html(	t.substring(0, start) +
									select_start +
									t.substring(start, end) +
									select_end +
									t.substring(end)
								);
	},
	get_sel: function()
	{
		if(!this._selected)
		{
			return ""
		}
		return this.seq.substring(this._select_start, this._select_end);
	},
	get_rev: function() {return this._reverse(this.get_sel());},
	get_cmp: function() {return this._complement(this.get_sel());},
	get_rev_cmp: function() {return this._reverse_complement(this.get_sel());},
	_update_toolbar: function()
	{
		if(this._selected)
		{
			$('.select_reqd').each( function(){
				$(this).button('enable');
			});
		}
		else
		{
			$('.select_reqd').each( function(){
				$(this).button('disable');
			});
		}
	},
//MOUSE FUNCTIONS
	_on_mouse_down: function(event)
	{
		this._mouse_down_stamp = event.timeStamp;
		this._mouse_down_pos = this._get_mouse_pos(event);
		this._active_row = this._get_mouse_row(event);
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
		//if(this._selected) this._clear_selection();
		this._select(this._mouse_down_pos, pos);
	},
	_window_mouse_up: function()
	{
		this._mouse_down = false;
	},
	_get_mouse_row: function(event)
	{
		return parseInt(event.currentTarget.id.split('-')[0]);
	},
	_get_mouse_pos: function(event) //given a mouse event, return the offset of the clicked character.
	{
		var largenum = parseInt(event.currentTarget.id.split('-')[0]);
		
		var smallnum = Math.floor((event.pageX - $(event.currentTarget).offset().left) / this.char_width);
		smallnum = toUnPadded(smallnum);
		
		return largenum + smallnum;
	},
	_on_scroll: function(event)
	{
		
		if($(window).scrollTop() > $('#seq_toolbar_rest').offset().top)
			$('#seq_toolbar').addClass('fixed').css('width', $('#seq_toolbar_rest').width());
		else
			$('#seq_toolbar').removeClass('fixed');
	},
	_get_char_width: function()
	{
		$r = $('#0-fwd');
		
		if($r.length == 0) this.char_width = char_width;
		else
		{
			this.char_width = $r.width() / toPadded(this.rowlength);
		}
		return this.char_width;
	},
//SEQUENCE FUNCTIONS
	_complement: function(string){
		var ret = "";
		for(i in string)
		{
			ret = ret + libFrag.alphabet[string[i]];
		}
		return ret;
	},
	_reverse: function(string){
		var ret = "";
		for(var i = string.length - 1; i >= 0; i = i-1)
		{
			ret = ret + string[i];
		}
		return ret;
	},
	_reverse_complement: function(string){
		var ret = "";
		for(var i = string.length - 1; i >= 0; i = i-1)
		{
			ret = ret + libFrag.alphabet[string[i]];
		}
		return ret;
	},
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
	if((f_e < (r_s+1)) || (f_s > r_e))
	{
		return "";
	}
	
	var left = f_s - r_s;
	if(left < 0) left = 0;
	var right = r_e - f_e;
	if(right < 0) right = 0;
	
	var l = r_e - r_s;
	
	var start = toPadded(left);
	var end = toPadded(l - right);
	l = toPadded(l);
	
	var r = "";
	for(var i = 0; i < l; i = i + 1)
	{
		if(i == start)
			r = r + '<span id="' + f_id + '-'+r_s+'-feat" class="feat-hl feat-type-' + feature.type.toLowerCase() + ' feat-id-' + f_id + '">';
		else if(i == end)
			r = r + '</span>';
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
	'<table style="word-wrap:break-word;width:100%;">' +
	'	<tr><td style="min-width:7.5em;"><h4>Type:</h4></td><td style="width:70%;"><p>' + feature.type + '</p></td></tr>' +
	'	<tr><td><h4>Location:</h4></td><td><p>' + (feature.start + 1) + '-' + (feature.end + 1) + '</p></td></tr>' + //NB. Added one to switcg to biologist-friendly 1-offset
	'	<tr><td><h4>Length:</h4></td><td><p>' + (feature.end - feature.start) + '</p></td></tr>' +
	'	<tr><td><h4>Strand:</h4></td><td><p>' + strand + '</p></td></tr>';
	
	for(i in feature.qualifiers)
	{
		ret = ret + '<tr><td><h4>' + feature.qualifiers[i].name + ':</h4></td><td><p class="feat-qual">' + feature.qualifiers[i].data + '</p></td></tr>';
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
