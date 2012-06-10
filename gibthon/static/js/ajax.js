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
	this.get = function(url, data, success_fn, fail_fn, async)
	{
	}

	//make an AJAX post request, calling success_fn or fail_fn with the returned
	//data
	this.post = function(url, data, success_fn, fail_fn, async)
	{
	}

	//Make a streaming AJAX request, i.e. call update_fn when each chunk of data
	//arrives and success_fn when all the data has arrived
	this.stream = function(url, data, success_fn, update_fn, fail_fn, method)
	{
	}

	/*
	 * Private functions
	 */

	var ajax_request(method, url, data, success_fn, fail_fn, async, update_fn)
	{
	}

}


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
