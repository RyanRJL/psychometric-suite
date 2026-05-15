/* Service worker for Psychometric Assistant.
   Strategy: cache-first with network fallback. On first visit the precache
   list is fetched. Subsequent visits load instantly from cache and update in
   the background. Any other requests (e.g. cache-busted assets) are cached
   on first fetch.

   Bump CACHE_VERSION below whenever you ship changes you want to force-
   refresh — old caches will be deleted on activation. */

const CACHE_VERSION = 'psyassist-v1';

/* Bare URLs (no ?v=...) — the SW also caches versioned variants on demand
   via the fetch handler below, so this list is just for first-paint speed. */
const PRECACHE_URLS = [
  './',
  './index.html',
  './styles.css',
  './design-system.css',
  './data.js',
  './app.js',
  './design-system.js',
  './app-effectsize-page.js',
  './manifest.json',
  './icon.svg',
  './icon-maskable.svg'
];

self.addEventListener('install', event => {
  event.waitUntil(
    caches.open(CACHE_VERSION)
      .then(cache => cache.addAll(PRECACHE_URLS).catch(err => {
        /* If any one URL 404s, don't abort the whole install — log and
           continue. The runtime fetch handler will pick up the missing
           file when it's actually requested. */
        console.warn('[sw] partial precache:', err);
      }))
      .then(() => self.skipWaiting())
  );
});

self.addEventListener('activate', event => {
  event.waitUntil(
    caches.keys()
      .then(keys => Promise.all(
        keys.filter(k => k !== CACHE_VERSION).map(k => caches.delete(k))
      ))
      .then(() => self.clients.claim())
  );
});

self.addEventListener('fetch', event => {
  /* Only GET requests are cacheable. Don't intercept anything else. */
  if (event.request.method !== 'GET') return;

  const url = new URL(event.request.url);
  /* Same-origin only — leave Google Fonts and any other CDN to their own
     HTTP caching. */
  if (url.origin !== location.origin) return;

  event.respondWith(
    caches.match(event.request, { ignoreSearch: true }).then(cached => {
      if (cached) return cached;
      return fetch(event.request).then(response => {
        if (response && response.ok){
          const clone = response.clone();
          caches.open(CACHE_VERSION).then(cache => cache.put(event.request, clone));
        }
        return response;
      }).catch(() => {
        /* Offline and not cached — for navigations, fall back to the
           shell so the SPA still loads. */
        if (event.request.mode === 'navigate'){
          return caches.match('./index.html');
        }
        return new Response('Offline', { status: 503, statusText: 'Offline' });
      });
    })
  );
});
