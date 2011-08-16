
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
				self.$origin.text("From: " + data[1].origin);
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
		this.$dd.slideDown(100);
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
		this.$dd.slideUp(100);
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
	_make_reference: function(ref):
	{
		
	},
});

})( jQuery );	
