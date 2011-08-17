
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
'	<div id="seq_wrap">' +
'	</div>' +
'</div>';


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
		this.$seq = this.$el.find('#seq_wrap');
		this.len = 0;
		this.prog = 0;
		this._get_alphabet();	
	},
	_get_seq: function(){
		var self = this;
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'seq', 'offset': this.prog}, function(data) {
			if(data[0] != 0)
			{
				console.log(data[1] + " while getting seq");
				return;
			}
			var olen = data[1].length;
			var oprog = self.prog;
			self.prog = self.prog + olen;
			
			if((self.prog) < self.len)
			{
				//we need to issue another request
				self._get_seq();			
			}
			else
			{
				self.$loader.slideUp(100);
			}
			
			added = 0;
			while(added < olen)
			{
				self.$seq.append(self._make_block(data[1], oprog + added));
				data[1] = data[1].slice(10);
				added = added + 10;
			}
			
			self.$prog.text(self.prog);
			self.$bar.progressbar('value', parseInt((100 * self.prog) / self.len));
		});
	},
	_get_alphabet: function(){
		var self = this;
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'alpha',}, function(data) {
			if(data[0] != 0)
			{
				console.log(data[1] + " while getting alphabet");
				return;
			}
			self.alphabet = data[1];
			self._get_len();
		});
	},
	_get_len: function(){
		var self = this;
		$.getJSON("/fragment/get/" + this.options.id + "/", {'value': 'len',}, function(data) {
			if(data[0] != 0)
			{
				console.log(data[1] + " while getting length");
				return;
			}
			self.len = data[1];
			self.$len.text(data[1]);
			if(self.len > 2000) //if we should load progressively 
			{
				self.$loader.slideDown(100);
			}
			self._get_seq();
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
	_make_block: function(sequence, start){
		seq = sequence.substr(0,10);
		cseq = this._complement(seq);
		return '' +
		'<div id="fe-content-' + start + '" class="fe-content-div">' +
		'	<div id="fe-feattop-' + start + '" class="fe-features">' +
		'	</div>' +
		'	<div style="position:relative;">'  +
		'		<div id="fe-htop-' + start + '" class="fe-hidden">' + seq + '</div>' +
		'		<div id="fe-hbottom-' + start + '" class="fe-hidden">' + cseq + '</div>' +
		'		<div id="fe-number-' + start + '" class="fe-label fe-number">' + start + '</div>' + 
		'		<div id="fe-marker-' + start + '" class="fe-label fe-marker">|</div>' +
		'		<div id="fe-top-' + start + '" class="fe-sequence fe-top">' + seq + '</div>' +
		'		<div id="fe-bottom-' + start + '" class="fe-sequence fe-bottom">' + cseq + '</div>' +
		'	</div>' +
		'	<div id="fe-featbottom-' + start + '" class="fe-features"></div>' + 
		'</div>';
	},
});

})( jQuery );	
