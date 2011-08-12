(function( $ ){
	$.fn.slideView = function () {
		return this.each(function() {
			var stop = false;
			var $this = $(this);
			var prefix = $this[0].id.replace('_wrapper','');
			$('#'+prefix+'_first')
				.button({
					icons:{primary:'ui-icon-seek-first'},
					disabled:true
				})
				.click(function(event, ui){
					if (stop){ return false;}
					var vis = $('div.'+prefix+':visible');
					var first = $('div.'+prefix+':first');
					stop = true;
					vis.hide('slide',{direction:"right"}, function() {$(first).show('slide',{direction:"left"}, function() {stop=false});});
					$(this).button('disable');
					$('#'+prefix+'_prev').button('disable');
					$('#'+prefix+'_next').button('enable');
					$('#'+prefix+'_last').button('enable');
				});
			$('#'+prefix+'_prev')
				.button({
					icons:{primary:'ui-icon-seek-prev'},
					disabled:true
				})
				.click(function(event, ui){
					if (stop){ return false;}
					var vis = $('div.'+prefix+':visible');
					var pre = $(vis[0].previousElementSibling);
					stop = true;
					vis.hide('slide',{direction:"right"}, function() {$(pre).show('slide',{direction:"left"}, function() {stop=false});});
					if (pre[0] == $('div.'+prefix+':first')[0]){
						$(this).button('disable');
						$('#'+prefix+'_first').button('disable');
					}
					$('#'+prefix+'_next').button('enable');
					$('#'+prefix+'_last').button('enable');
				});
			$('#'+prefix+'_last')
				.button({
					icons:{secondary:'ui-icon-seek-end'},
				})
				.click(function(event, ui){
					if (stop){ return false;}
					var vis = $('div.'+prefix+':visible');
					var last = $('div.'+prefix+':last');
					stop = true;
					vis.hide('slide',{direction:"left"}, function() {$(last).show('slide',{direction:"right"}, function() {stop=false});});
					$(this).button('disable');
					$('#'+prefix+'_next').button('disable');
					$('#'+prefix+'_prev').button('enable');
					$('#'+prefix+'_first').button('enable');
				});
			$('#'+prefix+'_next')
				.button({
					icons:{secondary:'ui-icon-seek-next'}
				})
				.click(function(event, ui){
					if (stop){ return false;}
					var vis = $('div.'+prefix+':visible');
					var nex = $(vis[0].nextElementSibling);
					stop = true;
					vis.hide('slide',{direction:"left"}, function() {$(nex).show('slide',{direction:"right"}, function() {stop=false});});
					if (nex[0] == $('div.'+prefix+':last')[0]){
						$(this).button('disable');
						$('#'+prefix+'_last').button('disable');

					}
					$('#'+prefix+'_prev').button('enable');
					$('#'+prefix+'_first').button('enable');
				});
				$('div.'+prefix+':first').show();
		});
	};
})( jQuery );