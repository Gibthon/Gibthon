

(function( $, undefined ) {

$.widget("ui.fragmentDrag", {
	options: {
		color: 'red',
	},
	_init: function()
	{
		var self = this;
		this.$el = $(this.element[0]);
		this.$fragments = new Array();
		$.getJSON('/fragment/list/', {}, function(data) {
			if(data[0] != 0) 
			{
				console.log('ERROR!');
				return;
			}
			for(var i in data[1].fragments)
			{
				var f = data[1].fragments[i];
				self._addFragment(f);
			}
			self._display();
		});
	},
	used: function(fid) //set the fragment fid to in use
	{
		this.$el.find('#frag_' + fid)
			.addClass('ui-state-disabled')
			.draggable('disable');
	},
	unused: function(fid) //set the fragment fid to in use
	{
		this.$el.find('#frag_' + fid)
			.removeClass('ui-state-disabled')
			.draggable('enable');
	},
	_addFragment: function(f)
	{
		var desc = f.desc.replace(/'/g, '"');
		console.log(desc);
		var $f = $("<div id='frag_" + f.id + "' class='fragment ui-widget ui-state-default ui-corner-all' title='" + desc + "' ><span class='fragment_name'>" + f.name + "</span></div>");
		this.$fragments.push($f.draggable({ 
			revert: true,
			containment: 'document',
			start: function() {$f.addClass('ui-state-highlight');},
			stop: function() {$f.removeClass('ui-state-highlight');},
		}));
	},
	_display: function()
	{
		var self = this;
		this.$el.fadeOut('fast', function() {
			$(this).html('');
			for(var i in self.$fragments)
			{
				$(this).append(self.$fragments[i]);
			}
			$(this).fadeIn('fast');
		});
	},
});

})(jQuery);
