/*
 * ui.loader
 * 
 * 	Performs a list of AJAX actions while displaying a progress bar
 *  Displays the error message if one accurs.
 *  Actions which return values (other than OK or ERROR) are not supported
 * 
 * usage: $(SELECTOR).loader({commands: [{	'desc': 'Command Description',
 * 											'url': '/url/to/call/',
 * 											'data': {	'the': 'data to post',
 * 														'to': 'the url'}, 
 * 											}, MORE COMMANDS HERE];
 * 
 * then call .loader('start') to begin, or set option autoStart to true.
 */

var loader_html = '' + 
'<div style="text-align:center;margin-top:1.5em;">' +
'	<img src="{{ STATIC_URL }}/images/spinner.gif" alt="Loading" />' +
'	<div id="command_list_div">' +
'		<ul id="command_list">' +
'		</ul>' +
'   </div>' +
'	<div id="final_msg" style="height:1.1em;margin:.2em 0;"></div>' +
'	<div id="loader_progress"></div>' +
'	<div id="buttons" style="display:none;"><button id="close_btn"></button></div>' +
'</div>';

(function( $, undefined ) {

$.widget("ui.loader", {
	options: {
		commands: [], //a list of commands, in the form [{'desc': DESCRIPTION, 'url': URL, 'data': DATA}, etc]
		autostart: false, //should the actions start immediately?
	},
	_create: function()
	{
		if(this.options.commands == [])
			this._reload_page();
		
		var self = this;
		this.$el = $(this.element[0]);
		this.$progress = this.$el.find('#loader_progress').progressbar({value: 0,});
		this._cmd = this.options.commands;
		this._len = this.options.commands.length;
		this._next = 0; //the offset of the function to be performed next
		this._errors = 0;
		var $cmdList = this.$el.find('#command_list');
		for(i in this._cmd)
		{
			$cmdList.append('<li id="command-' + i + '" class="command_li"><p>' + this._cmd[i].desc + '<p></li>');
		}
		this._started = false;
		if(this.options.autoStart)
		{
			this._started = true;
			this._next_action();
		}
	},
	start: function() //start performing the commands
	{
		if(!this._started)
		{
			this._started = true;
			this._next_action();
		}
	},
	progress: function() //return the progress
	{
		return this._next;
	},
	_next_action: function()
	{
		//perform the next action
		if(this._next >= this._len)
			this._on_complete();
		var self = this;
		var cmd = this._cmd[this._next];
		$.getJSON( cmd.url, cmd.data, function(data){
			if(data[0] != 0)
			{
				//there was an error! Show the error text
				$('#command-' + self._next)
					.removeClass('currentCommand')
					.addClass('errorCommand')
					.append('<div class="errorDetail"><p>' + data[1] + '</p></div>');
				this._errors = this._errors + 1;
			}
			else //everything OK
			{
				$('#command-' + self._next)
					.removeClass('currentCommand')
					.addClass('okCommand');
			}
			//move on to the next command
			self._next = self._next + 1;
			$('#command-' + self.next).addClass('currentCommand');
			self._update_progressbar();
			self._next_action();
		});
	},
	_on_complete: function()
	{
		if(this._errors > 0)
		{
			var e = ' error';
			if(this._errors > 1) e = e + 's';
			this.$el.find('#final_msg').text('There were ' + this._errors + e + '.');
			$('#close_btn').button({
				label: 'Close',
				icons: {primary: ''},
			}).click( function(event) {
				self._reload_page();
			});
			$('#buttons').slideDown(100);
			return;
		}
		self._reload_page(3000); //wait 3 secs, then reload
		
	},
	_update_progressbar: fucntion()
	{
		this.$progress.progressbar('value', this._next / this._len);
	},	
	_reload_page: function(delay)
	{
		//some JS to reaload the page after a delay
		if(delay == undefined) delay = 0;
		console.log("TODO: FIND OUT HOW TO RELOAD THE PAGE AFTER A DELAY...");
		
	}
});

})(jQuery);
