var inbox_refresh = function() {
	$.get('/user/inbox/fetch', function (data) {
		$('span#inbox_not_added').html(data[2]);
		$('span#inbox_unread').html(data[1]);
		$('span#inbox_count').html(data[3]);
		if (data[0] > 0) {
			alert('You have '+data[0]+' new message(s)!');
		}
	}, 'json');
};

$(document).ready(function() {
	inbox_refresh();
			
	if ($('#user_tab_profile').length == 1) {
		$('#user_tab_profile').button({
			icons:{primary:'ui-icon-person',
				   secondary:'ui-icon-triangle-1-s'}
		});
		$('#user_tab_inbox').button({
			icons:{primary:'ui-icon-mail-closed'}
		});
	} else {
		$('#user_tab_login').button({
			icons:{primary:'ui-icon-tag'}
		})
		.click(function (event) {
			event.preventDefault();
			window.location.href="/user/login?next="+window.location.pathname;
		});
		$('#user_tab_register').button({
			icons:{primary:'ui-icon-key'}
		});
	}
	$('#user_tab_home').button({
		icons:{primary:'ui-icon-home',
			   secondary:'ui-icon-triangle-1-s'}
	});
	$('#user_tab_help').button({
		icons:{primary:'ui-icon-help',
			   secondary:'ui-icon-triangle-1-s'}
	});
	$('#user_tab_tools').button({
		icons:{primary:'ui-icon-wrench',
			   secondary:'ui-icon-triangle-1-s'}
	});
	
	$('.user_tab').addClass('ui-corner-top').removeClass('ui-corner-all');

	$('span.drop').mouseleave(function () {
		$('ul', this).slideUp("fast") }
	);
	
	 $('span.drop a.user_tab').mouseenter(function () {$('~ ul', this).slideDown('fast')});
	
	$('span.drop ul li a').button();
	
	$('input[id^="id_captcha"]').addClass('captcha');
});
  var _gaq = _gaq || [];
  _gaq.push(['_setAccount', 'UA-18410656-4']);
  _gaq.push(['_trackPageview']);

  (function() {
    var ga = document.createElement('script'); ga.type = 'text/javascript'; ga.async = true;
    ga.src = ('https:' == document.location.protocol ? 'https://ssl' : 'http://www') + '.google-analytics.com/ga.js';
    var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(ga, s);
  })();
