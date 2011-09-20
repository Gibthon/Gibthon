/*
 * jQuery UI Prompt
 *
 * Displays an alert using the jQuery UI theme
 * 
 * Written by Bill Collins for Gibthon
 */

function makeHttpObject() {
  try {return new XMLHttpRequest();}
  catch (error) {}
  try {return new ActiveXObject("Msxml2.XMLHTTP");}
  catch (error) {}
  try {return new ActiveXObject("Microsoft.XMLHTTP");}
  catch (error) {}

  throw new Error("Could not create HTTP request object.");
}

 function getCookie(name) {
	var cookieValue = null;
	if (document.cookie && document.cookie != '') {
		var cookies = document.cookie.split(';');
		for (var i = 0; i < cookies.length; i++) {
			var cookie = jQuery.trim(cookies[i]);
			// Does this cookie string begin with the name we want?
			if (cookie.substring(0, name.length + 1) == (name + '=')) {
				cookieValue = decodeURIComponent(cookie.substring(name.length + 1));
				break;
			}
		}
	}
	return cookieValue;
}

(function ($, undefined) {

$.widget( "ui.prompt" , {
	options: {
		title: "Alert",
		message: "Body",
		modal: true,
		type: 'alert',
		confirm: {
			text: "OK",
			click: function (promptText) { return true },
			icon: "ui-icon-check"
		},
		cancel: {
			text: "Cancel",
			click: function () { return false },
			icon: "ui-icon-close"
		},
		promptText: '',
		target: { 
			location :'',
			target: '',
			type: 'GET',
			data: '',
			callback: '',
		},
		buttons: true
	},
	_create: function() {
		var buttons;
		var message = this.options.message;
		var self = this;
		var confirmButton = { "text":this.options.confirm.text, "click": function () { self._confirm() } };
		var cancelButton = { "text":this.options.cancel.text, "click": function () { self._cancel() } };
		console.log(this.options.confirm);
		switch(this.options.type) {
			case 'alert':
				buttons = [ confirmButton ];
				break;
			case 'confirm':
				buttons = [ cancelButton, confirmButton ];
				break;
			case 'prompt':
				buttons = [ cancelButton, confirmButton ];
				message = message + '<br /><br /><input type="text" style="width:90%;" value="'+this.options.promptText+'" id="promptText" />';
				break;
			case 'wait':
				message = message + '<br /><br /><div style="text-align:center;"><img src="'+STATIC_URL+'/images/spinner.gif" /></div>';
				break;
			case 'progress':
				if (this.options.target.location == '') throw("No target specified")
				buttons = [ cancelButton, confirmButton ];
				message = message + '<br /><br /><div id="progress" style="height:22px;"></div>'
				break;
			default:
				throw("Invalid dialog type");
		}
		if (!(this.options.buttons)) {
			buttons = '';
		}
		this.element
			.html(message)
			.dialog({
				title: this.options.title,
				resizable: false,
				closeOnEscape: false,
				draggable: false,
				modal: this.options.modal,
				buttons: buttons,
				autoOpen:false,
				open: function(event, ui) { $(".ui-dialog-titlebar-close").hide() },
				show: "fade",
				hide: "fade"
			});
		
		switch($('~ div button', this.element).size()) {
			case 0:
				break;
			case 1:
				if (this.options.confirm.icon != "") {
					$('~ div button:first', this.element).button("option", "icons", {primary:this.options.confirm.icon});
				}
				break;
			case 2:
				if (this.options.confirm.icon != "") {
					$('~ div button:nth-child(2)', this.element).button("option", "icons", {primary:this.options.confirm.icon});
				}
				if (this.options.cancel.icon != "") {
					$('~ div button:first-child', this.element).button("option", "icons", {primary:this.options.cancel.icon});
				}
				break;
			default:
		}
		this.element.dialog('open');
		if (this.options.type == "progress") {
			$('~ div button', this.element).button("disable");
			$('#progress', this.element).progressbar({ value:0 });
			var h = makeHttpObject();
			h.open(this.options.target.type, this.options.target.location, true)
			h.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
			h.onreadystatechange = function (){
				p = parseInt(this.responseText.split(':').pop());
				$('#progress', self.element).progressbar('value', p);
				if (p==100){
					$('~ div button', self.element).button("enable");
					if( self.options.target.callback != '' ) self.options.target.callback();
				}
			}
			h.send(this.options.target.data);
		}
		if (this.options.type == "wait") {
			$(this.options.target.target).load(this.options.target.location, function () { self.options.target.callback(); self.close(); });
		}
		
	},
	
	_text: function () {
		return $('input', this.element).val();
	},
	
	_confirm: function () {
		if (!(this.options.type == 'prompt' && this._text() == '')) {
			this.options.confirm.click(this._text());
			this.close();
		} else {
			this._cancel();
		}
	},
	
	_cancel: function () {
		this.options.cancel.click();
		this.close();
	},

	
	close: function () {
		this.element.dialog("close");
		this._destroy();
	},
	
	_destroy: function() {
		this.element.html('');
		this.element.dialog('destroy');
		$.Widget.prototype.destroy.apply(this, arguments);
	}
});
})(jQuery);