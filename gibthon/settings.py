# Django settings for gibthon project.

from localsettings import localsettings

# must contain the following values:
#	email_password
#	db_password
#	db_name
#	db_user
#	key (django specific key)
#	unafold_wd  (path to temp directory for unafold files)
#	hybrid_ss_min_path (path to hyrbid-ss-min)
#	hybrid_min_path (path to hybrid-min)
#	boxplot_ng_path (path to boxplot_ng)
#	template_dirs (tuple of template directories, absolute paths)
#	static_dirs (tuple of static directories, absolute paths)
#	static_root (static file collection location) 
#	media_root (media file location)
#	media_url (media url)
#	static_url (static url)

DEBUG = localsettings.debug
TEMPLATE_DEBUG = DEBUG

ADMINS = (
    ('Bill Collins', 'bill@gibthon.org'),
		('Haydn King', 'haydn@gibthon.org'),
)

MANAGERS = ADMINS

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql', # Add 'postgresql_psycopg2', 'postgresql', 'mysql', 'sqlite3' or 'oracle'.
        'NAME': localsettings.db_name,                      # Or path to database file if using sqlite3.
        'USER': localsettings.db_user,                      # Not used with sqlite3.
        'PASSWORD': localsettings.db_password,                  # Not used with sqlite3.
        'HOST': '',                      # Set to empty string for localhost. Not used with sqlite3.
        'PORT': '',                      # Set to empty string for default. Not used with sqlite3.
    }
}

# Local time zone for this installation. Choices can be found here:
# http://en.wikipedia.org/wiki/List_of_tz_zones_by_name
# although not all choices may be available on all operating systems.
# On Unix systems, a value of None will cause Django to use the same
# timezone as the operating system.
# If running in a Windows environment this must be set to the same as your
# system time zone.
TIME_ZONE = 'Europe/London'

# Language code for this installation. All choices can be found here:
# http://www.i18nguy.com/unicode/language-identifiers.html
LANGUAGE_CODE = 'en-gb'

SITE_ID = 1
DEFAULT_CHARSET = 'utf-8'

# If you set this to False, Django will make some optimizations so as not
# to load the internationalization machinery.
USE_I18N = False

# If you set this to False, Django will not format dates, numbers and
# calendars according to the current locale
USE_L10N = True

# Absolute filesystem path to the directory that will hold user-uploaded files.
# Example: "/home/media/media.lawrence.com/"
MEDIA_ROOT = localsettings.media_root

# URL that handles the media served from MEDIA_ROOT. Make sure to use a
# trailing slash if there is a path component (optional in other cases).
# Examples: "http://media.lawrence.com", "http://example.com/media/"
MEDIA_URL = localsettings.media_url

STATIC_ROOT = localsettings.static_root

STATIC_URL = localsettings.static_url

STATICFILES_DIRS = localsettings.static_dirs


# URL prefix for admin media -- CSS, JavaScript and images. Make sure to use a
# trailing slash.
# Examples: "http://foo.com/media/", "/media/".
ADMIN_MEDIA_PREFIX = '/static/admin/'

# Make this unique, and don't share it with anybody.
SECRET_KEY = localsettings.key

# List of callables that know how to import templates from various sources.
TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
#     'django.template.loaders.eggs.Loader',
)

MIDDLEWARE_CLASSES = (
    'django.middleware.common.CommonMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
)

ROOT_URLCONF = 'gibthon.urls'

TEMPLATE_DIRS = localsettings.template_dirs


INSTALLED_APPS = (
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    #'django.contrib.sites',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'django.contrib.admin',
    'django.contrib.admindocs',
		'south',
		'fragment',
    'gibson',
    'molcal',
    'captcha',
    'partsearch',
    'ligcalc',
    'accounts',
    'pricalc',
    'digest',
)

LOGIN_URL = '/user/login'
LOGIN_REDIRECT_URL = '/user/profile/'

AUTHENTICATION_BACKENDS = (
	'gibthon.auth_backends.GibthonUserModelBackend',
)

CUSTOM_USER_MODEL = 'accounts.GibthonUser'

# captcha settings

CAPTCHA_CHALLENGE_FUNCT = 'captcha.helpers.math_challenge'

# email settings

#EMAIL_HOST = 'smtp.gibthon.org'
#EMAIL_PORT = 25
#EMAIL_HOST_USER = 'site@gibthon.org'
#EMAIL_HOST_PASSWORD = localsettings.email_password
AWS_ACCESS_KEY_ID = localsettings.ses_accesskey
AWS_SECRET_ACCESS_KEY = localsettings.ses_secretaccesskey
EMAIL_BACKEND = 'django_ses.SESBackend'
EMAIL_USE_TLS = False
EMAIL_SUBJECT_PREFIX = '[gibthon.org] '
SERVER_EMAIL = 'admin@gibthon.org'

# local settings

UNAFOLD_WD = localsettings.unafold_wd
HYBRID_SS_MIN_PATH = localsettings.hybrid_ss_min_path
HYBRID_MIN_PATH = localsettings.hybrid_min_path
BOXPLOT_NG_PATH = localsettings.boxplot_ng_path
