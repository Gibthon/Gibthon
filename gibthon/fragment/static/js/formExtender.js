/*
formExtender:
* expands a single form field into multiple form fields
* 
* initialHTML
<div id='extender'> // the form extender
	<div id='extender-copy' style='display:none'> 
		everything in here is copied
		<button class='extender-remove'></button>
	</div>
	<div id='extender-target'></div>
	
	<button class='extender-add'></button>
</div>
*/

(function( $, undefined ) {

$.widget("ui.formExtender", {
	options: {
		addInitial: true, //whether to add an initial element
		},
	_init: function() {
		var self = this;
		this.$el = $(this.element[0]);
		
	//intialise the elements	
		this.$copy = this.$el.find('#extender-copy').detach(); //remove from the dom - it'll only cause trouble
		this.$target = this.$el.find('#extender-target');
		
		this.$el.find('button.extender-add').button({
			label: 'More',
			icons: {primary: 'ui-icon-plusthick',},
		}).click( function(event) {
			self._add(event);
		});
	//intialise the data
		this.num_elements = 0;
		
		if(this.options.addInitial)
		{
			this._add();
		}
	},
	numElements: function() {return num_elements},
	_add: function(event){
		var self=this;
		var $clone = this.$copy.clone();
		
		$clone.find('button.extender-remove').button({
			text: false,
			icons: {primary: 'ui-icon-trash',},
		}).click( function(event) {
			console.log('remove');
			self._remove(event);
		});
		
		$clone.appendTo(this.$target)
			.slideDown('fast');
		this.num_elements = this.num_elements + 1;
	},
	_remove: function(event){
		if(this.num_elements <= 1)
			return;
		$(event.target).closest('#extender-copy').slideUp('fast', function() {
			$(this).remove();
		});
		this.num_elements = this.num_elements - 1;
	},
});

})(jQuery);
