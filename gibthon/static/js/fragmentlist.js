/*
*	Javascript to power the list selection
* 
*  bind to fragmentList.selectChanged to get updated selection info
*  getList returns a list of all the selected fragments
* 
* 
*  MUST BE APPLIED BEFORE DATATABLES
* 
*/

(function( $, undefined ) {

$.widget("ui.fragmentList", {
	options: {
	},
	_init: function() {
		var $element = $(this.element[0]);
		var self = this;
		//bind to all the onclicks
		this.fragments = $element.find('.selected-check')
			.click(function(event) {self._update(event);});
			
		this.num_fragments = this.fragments.length;
		this.num_selected = 0;

		this.all_selected = $element.find('.all-selected')
			.click(function(event) {self._select_all(event);});
		
	},
	getNumSelected: function(){
		return this.num_selected;
	},
	_update: function(event) {
		//update the number of ones selected	
		this.num_selected = this.fragments.filter(':checked').length;		
		if(this.num_selected == this.num_fragments)
			this.all_selected.prop('checked', true);
		else
			this.all_selected.prop('checked', false);
		
		this._trigger('selectChanged', event, {'num': this.num_fragments,
															'selected': this.num_selected});
	},
	_select_all: function(event) {
		//select them all, or not
		if(this.all_selected.prop('checked'))
		{
			this.fragments.prop('checked', true);
		}
		else
		{
			this.fragments.prop('checked', false);
		}
		this.num_selected = this.fragments.filter(':checked').length;		
		this._trigger('selectChanged', event, {'num': this.num_fragments,
															'selected': this.num_selected});
	},
	getList: function() {
		var ret = new Array();
		this.fragments.filter(':checked').each(function(index, element) {
			ret.push($(element).prop('value'));
		});
		return ret;
	}
});

})( jQuery );	
