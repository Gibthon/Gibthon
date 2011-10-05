
(function( $, undefined ) {

$.widget("ui.designer", {
	options: {
		width: 10,
		height: 10,
	},
	_init: function()
	{
		var self = this;
		this.canvas = this.element[0];
		$(this.canvas).css({
			'width': this.options.width,
			'height': this.options.height,
		});
		
	}
	
});

})(jQuery);
