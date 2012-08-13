/* 
 * New AJAX API
 *   -- wraps around jQuery and the Gibthon standard for AJAX error reporting
 *      (i.e. JSONResponse)
 *
 */

var AJAX = new function()
{
	/*
	 * Public API
	 *
	 */

	//Make an AJAX get request, and return to success_fn on success, or fail_fn on
	//failure.
	this.get = function(args)
	{
		args.type = 'GET';
		ajax_request(args);
	}

	//make an AJAX post request, calling success_fn or fail_fn with the returned
	//data
	this.post = function(args)
	{
		args.type = 'POST';
		ajax_request(args);
	}

	//Make a streaming AJAX request, i.e. call update_fn when each chunk of data
	//arrives and success_fn when all the data has arrived
	this.stream = function(args, update_fn)
	{
		//Calling without an update_fn is pointless -- you should use get or post
		//above
		if(update_fn==undefined)
		{
			console.error('AJAX.stream: called without update function');
			//treat it as a normal request
			ajax_request(args);
			return;
		}

		var xhr = makeHttpRequest();

		xhr.onreadystatechange = function()
		{
			//If we're done
			if(this.readyState==4)
				{
					//and there are no errors
					if(this.status==200)
						{
							if(args.success!=undefined)
								args.success(this.responseText);
						}
					else
						{
							if(args.error!=undefined)
								args.error(this.statusText)
						}
				}
				//if more data arrived
				if(this.readyState==3)
					update_fn(this.responseText);
		}

		//open and send
		if(args.type==undefined) args.type = 'POST';
		if(args.url==undefined) args.url = '';
		xhr.open(args.type, args.url);
		xhr.send(args.data);
	}

	/*
	 * Private functions
	 */

	var ajax_request = function(args)
	{
		//replace the success function with my own
		var s_fn = args.success;
		args.success = function(data, textStatus, jqXHR)
		{
			//check for error 
			if(data[0] == AJAX_ERROR)
				{
					//call the error callback, if it exists
					if($.isFunction(args.error))
						args.error(jqXHR, data[1], 'ServerError');
					else
						console.error('AJAX Error: ' + data[1]);
				}
				else
					if($.isFunction(s_fn)) s_fn(data[1]);
		}

		//stringify the data correctly if it hasn't been already
		//if(args.data!=undefined)
			//args.data = JSON.stringify(args.data);

		$.ajax(args);
	}

	var makeHttpRequest = function() 
	{
		try 
		{
			return new XMLHttpRequest();
		}
		catch (error) {}
		
		try 
		{
			return new ActiveXObject("Msxml2.XMLHTTP");
		}
		catch (error) {}
		
		try 
		{
			return new ActiveXObject("Microsoft.XMLHTTP");
		}
		catch (error) {}

		throw new Error("Could not create HTTP request object.");
	}

	/*
	 * Private Variables
	 */

	//value returned when there's an error
	var AJAX_ERROR = -1;
};


/*
 * Legacy AJAX code
 *
 */

$(document).ajaxSend(function(event, xhr, settings) {
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
	function sameOrigin(url) {
		// url could be relative or scheme relative or absolute
		var host = document.location.host; // host + port
		var protocol = document.location.protocol;
		var sr_origin = '//' + host;
		var origin = protocol + sr_origin;
		// Allow absolute or scheme relative URLs to same origin
		return (url == origin || url.slice(0, origin.length + 1) == origin + '/') ||
			(url == sr_origin || url.slice(0, sr_origin.length + 1) == sr_origin + '/') ||
			// or any other URL that isn't scheme relative or absolute i.e relative.
			!(/^(\/\/|http:|https:).*/.test(url));
	}
	function safeMethod(method) {
		return (/^(GET|HEAD|OPTIONS|TRACE)$/.test(method));
	}

	if (!safeMethod(settings.type) && sameOrigin(settings.url)) {
		xhr.setRequestHeader("X-CSRFToken", getCookie('csrftoken'));
	}
});
