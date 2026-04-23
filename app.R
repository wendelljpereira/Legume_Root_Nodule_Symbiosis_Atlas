# Legume Root Nodule Symbiosis Atlas Shiny entrypoint.
# Reusable helpers and tab builders live in R/ so the app can be tested and
# maintained without keeping every layer in one monolithic file.
source(file.path("R", "atlas_core.R"), local = TRUE)
source(file.path("R", "atlas_ui_tabs.R"), local = TRUE)

ui <- fluidPage(
    tags$head(
        tags$link(rel = "preconnect", href = "https://fonts.googleapis.com"),
        tags$link(rel = "preconnect", href = "https://fonts.gstatic.com", crossorigin = "anonymous"),
        tags$link(
            rel = "stylesheet",
            href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=IBM+Plex+Sans:wght@500;600;700&family=IBM+Plex+Mono:wght@400;500&display=swap"
        ),
        tags$link(rel = "icon", type = "image/svg+xml", href = "favicon.svg"),
        tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
        tags$script(HTML(
            "(function() {
                function atlasVisiblePermalinkPanel() {
                    var anchors = Array.prototype.slice.call(document.querySelectorAll('.permalink-panel'));
                    anchors = anchors.filter(function(el) {
                        return el && el.offsetParent !== null;
                    });
                    if (!anchors.length) return null;

                    var targetOffset = 132;
                    var best = null;
                    var bestScore = Infinity;

                    anchors.forEach(function(el) {
                        var rect = el.getBoundingClientRect();
                        var visible = rect.bottom > targetOffset && rect.top < window.innerHeight * 0.8;
                        var score = Math.abs(rect.top - targetOffset);
                        if (!visible) score += 10000;
                        if (score < bestScore) {
                            bestScore = score;
                            best = el;
                        }
                    });

                    return best ? best.getAttribute('data-permalink-panel') : null;
                }

                function atlasFindPermalinkPanel(panelId) {
                    var anchors = Array.prototype.slice.call(document.querySelectorAll('.permalink-panel'));
                    for (var i = 0; i < anchors.length; i += 1) {
                        if (anchors[i].getAttribute('data-permalink-panel') === panelId) {
                            return anchors[i];
                        }
                    }
                    return null;
                }

                function atlasShowToast(message, tone) {
                    var toast = document.getElementById('atlas-toast');
                    if (!toast) {
                        toast = document.createElement('div');
                        toast.id = 'atlas-toast';
                        toast.className = 'atlas-toast';
                        document.body.appendChild(toast);
                    }

                    toast.textContent = message;
                    toast.className = 'atlas-toast is-visible ' + (tone || 'success');

                    if (toast._hideTimer) {
                        window.clearTimeout(toast._hideTimer);
                    }

                    toast._hideTimer = window.setTimeout(function() {
                        toast.className = 'atlas-toast';
                    }, 2200);
                }

                function atlasAppBusy() {
                    return !!document.querySelector('.shiny-busy, .recalculating');
                }

                window.atlasBulkGeneOptionValue = function(query) {
                    return '__bulk__::' + encodeURIComponent(query || '');
                };

                window.atlasUpdateBulkGeneOption = function(selectize, query) {
                    if (!selectize) return;

                    Object.keys(selectize.options || {}).forEach(function(key) {
                        if (key.indexOf('__bulk__::') === 0) {
                            selectize.removeOption(key, true);
                        }
                    });

                    query = (query || '').trim();
                    if (!query.length) {
                        selectize.refreshOptions(false);
                        return;
                    }

                    var optionValue = window.atlasBulkGeneOptionValue(query);
                    selectize.addOption({
                        value: optionValue,
                        label: 'Select all matches for \\\"' + query + '\\\"',
                        tokens: query,
                        bulk_select: true
                    });
                    selectize.refreshOptions(false);

                    window.requestAnimationFrame(function() {
                        var optionNode = selectize.getOption(optionValue);
                        if (optionNode && optionNode.length) {
                            optionNode.addClass('selectize-bulk-option');
                            var parent = optionNode.parent();
                            if (parent && parent.length) {
                                parent.prepend(optionNode);
                            }
                        }
                    });
                };

                window.atlasScheduleBulkGeneOption = function(selectize, query) {
                    window.atlasUpdateBulkGeneOption(selectize, query);
                    [90, 220].forEach(function(delayMs) {
                        window.setTimeout(function() {
                            window.atlasUpdateBulkGeneOption(selectize, query);
                        }, delayMs);
                    });
                };

                window.atlasHandleBulkGeneOption = function(selectize, value) {
                    if (!value || value.indexOf('__bulk__::') !== 0) {
                        return false;
                    }

                    var query = '';
                    try {
                        query = decodeURIComponent(value.replace('__bulk__::', ''));
                    } catch (e) {
                        query = value.replace('__bulk__::', '');
                    }

                    window.setTimeout(function() {
                        selectize.removeItem(value, true);
                        selectize.setTextboxValue('');
                        selectize.close();
                        window.atlasUpdateBulkGeneOption(selectize, '');
                    }, 0);

                    if (window.Shiny && Shiny.setInputValue) {
                        Shiny.setInputValue('selected_genes_bulk_query', {
                            query: query,
                            nonce: Date.now()
                        }, { priority: 'event' });
                    }

                    return true;
                };

                function atlasCloseOpenGenePickers() {
                    var activeControls = Array.prototype.slice.call(document.querySelectorAll('.selectize-control.multi.focus'));
                    activeControls.forEach(function(control) {
                        var input = control.querySelector('input[id$=\"-selectized\"]');
                        if (!input || !input.selectize) return;

                        input.selectize.setTextboxValue('');
                        input.selectize.close();
                        input.selectize.refreshOptions(false);
                    });
                }

                document.addEventListener('mousedown', function(event) {
                    if (event.target && event.target.closest && event.target.closest('.apply-selection-btn')) {
                        atlasCloseOpenGenePickers();
                    }
                }, true);

                document.addEventListener('touchstart', function(event) {
                    if (event.target && event.target.closest && event.target.closest('.apply-selection-btn')) {
                        atlasCloseOpenGenePickers();
                    }
                }, true);

                function atlasSetButtonBusy(buttonId, busy, busyLabel) {
                    var button = document.getElementById(buttonId);
                    if (!button) return;

                    if (!button.dataset.atlasOriginalHtml) {
                        button.dataset.atlasOriginalHtml = button.innerHTML;
                    }

                    if (busy) {
                        button.disabled = true;
                        button.classList.add('is-working');
                        button.setAttribute('aria-busy', 'true');
                        button.innerHTML = '<span class=\"atlas-inline-spinner\" aria-hidden=\"true\"></span><span>' + (busyLabel || 'Working...') + '</span>';
                        return;
                    }

                    button.disabled = false;
                    button.classList.remove('is-working');
                    button.removeAttribute('aria-busy');

                    if (button.dataset.atlasOriginalHtml) {
                        button.innerHTML = button.dataset.atlasOriginalHtml;
                    }
                }

                function atlasReportVisiblePermalinkPanel(panelId) {
                    if (!(window.Shiny && Shiny.setInputValue)) {
                        return;
                    }

                    if (!panelId && !window.__atlasRestoringPanelId && atlasAppBusy()) {
                        return;
                    }

                    Shiny.setInputValue('visible_permalink_panel', {
                        panel: panelId || window.__atlasRestoringPanelId || atlasVisiblePermalinkPanel() || '',
                        nonce: Date.now()
                    }, { priority: 'event' });
                }

                function atlasPinboardState() {
                    if (!window.__atlasPinboardState) {
                        window.__atlasPinboardState = {
                            items: [],
                            maxItems: 4,
                            nextId: 1
                        };
                    }
                    return window.__atlasPinboardState;
                }

                function atlasEnsurePinboardUi() {
                    var root = document.getElementById('atlas-pinboard');
                    if (!root) return null;

                    return {
                        root: root,
                        empty: root.querySelector('.atlas-pinboard-empty'),
                        grid: root.querySelector('.atlas-pinboard-grid'),
                        clear: root.querySelector('.atlas-pinboard-clear')
                    };
                }

                function atlasActiveTabLabel() {
                    var active = document.querySelector('.nav-tabs li.active a');
                    return active ? (active.textContent || '').trim() : 'Current tab';
                }

                function atlasSerializeSvg(svgEl) {
                    var svgMarkup = new XMLSerializer().serializeToString(svgEl);
                    return 'data:image/svg+xml;charset=utf-8,' + encodeURIComponent(svgMarkup);
                }

                function atlasCapturePlotCard(card) {
                    var titleNode = card.querySelector('.plot-card-title');
                    var title = titleNode ? (titleNode.textContent || '').trim() : 'Pinned plot';
                    var context = atlasActiveTabLabel();

                    return new Promise(function(resolve, reject) {
                        var plotlyEl = card.querySelector('.js-plotly-plot');
                        if (plotlyEl && window.Plotly && window.Plotly.toImage) {
                            window.Plotly.toImage(plotlyEl, {
                                format: 'png',
                                width: Math.max(plotlyEl.clientWidth || 0, 880),
                                height: Math.max(plotlyEl.clientHeight || 0, 520)
                            }).then(function(dataUrl) {
                                resolve({ kind: 'img', src: dataUrl, title: title, context: context });
                            }).catch(function() {
                                reject(new Error('plotly-export-failed'));
                            });
                            return;
                        }

                        var imgEl = card.querySelector('img');
                        if (imgEl && (imgEl.currentSrc || imgEl.src)) {
                            resolve({ kind: 'img', src: imgEl.currentSrc || imgEl.src, title: title, context: context });
                            return;
                        }

                        var canvasEl = card.querySelector('canvas');
                        if (canvasEl && canvasEl.toDataURL) {
                            resolve({ kind: 'img', src: canvasEl.toDataURL('image/png'), title: title, context: context });
                            return;
                        }

                        var svgEl = card.querySelector('svg');
                        if (svgEl) {
                            resolve({ kind: 'img', src: atlasSerializeSvg(svgEl), title: title, context: context });
                            return;
                        }

                        reject(new Error('plot-not-ready'));
                    });
                }

                function atlasRenderPinboard() {
                    var ui = atlasEnsurePinboardUi();
                    if (!ui) return;

                    var state = atlasPinboardState();
                    ui.grid.innerHTML = '';

                    if (!state.items.length) {
                        ui.root.classList.remove('has-items');
                        return;
                    }

                    ui.root.classList.add('has-items');

                    state.items.forEach(function(item) {
                        var tile = document.createElement('div');
                        tile.className = 'atlas-pinboard-tile';

                        var head = document.createElement('div');
                        head.className = 'atlas-pinboard-tile-head';

                        var textWrap = document.createElement('div');
                        textWrap.className = 'atlas-pinboard-tile-meta';

                        var title = document.createElement('div');
                        title.className = 'atlas-pinboard-title';
                        title.textContent = item.title;

                        var context = document.createElement('div');
                        context.className = 'atlas-pinboard-context';
                        context.textContent = item.context;

                        var removeBtn = document.createElement('button');
                        removeBtn.type = 'button';
                        removeBtn.className = 'atlas-pinboard-remove';
                        removeBtn.setAttribute('data-pinboard-remove', item.id);
                        removeBtn.setAttribute('aria-label', 'Remove pinned plot');
                        removeBtn.textContent = 'Remove';

                        textWrap.appendChild(title);
                        textWrap.appendChild(context);
                        head.appendChild(textWrap);
                        head.appendChild(removeBtn);

                        var preview = document.createElement('div');
                        preview.className = 'atlas-pinboard-preview';

                        if (item.kind === 'img') {
                            var img = document.createElement('img');
                            img.src = item.src;
                            img.alt = item.title;
                            preview.appendChild(img);
                        }

                        tile.appendChild(head);
                        tile.appendChild(preview);
                        ui.grid.appendChild(tile);
                    });
                }

                function atlasPinPlotCard(card) {
                    if (!card) return;

                    atlasCapturePlotCard(card).then(function(snapshot) {
                        var state = atlasPinboardState();
                        snapshot.id = 'pin-' + state.nextId;
                        state.nextId += 1;

                        if (state.items.length >= state.maxItems) {
                            state.items.shift();
                            atlasShowToast('Scratchpad full. Replaced the oldest plot.', 'success');
                        }

                        state.items.push(snapshot);
                        atlasRenderPinboard();
                        atlasShowToast('Pinned plot', 'success');
                    }).catch(function() {
                        atlasShowToast('Wait for the plot to finish rendering before pinning it.', 'error');
                    });
                }

                function atlasDecoratePlotCards() {
                    var cards = Array.prototype.slice.call(document.querySelectorAll('.plot-card'));

                    cards.forEach(function(card) {
                        if (card.__atlasPinReady) return;

                        var title = card.querySelector(':scope > .plot-card-header .plot-card-title, :scope > .plot-card-title');
                        if (!title) return;

                        var header = card.querySelector(':scope > .plot-card-header');
                        if (!header) {
                            header = document.createElement('div');
                            header.className = 'plot-card-header';
                            title.parentNode.insertBefore(header, title);
                            header.appendChild(title);
                        }

                        var pinBtn = document.createElement('button');
                        pinBtn.type = 'button';
                        pinBtn.className = 'plot-card-pin-btn plot-download-btn';
                        pinBtn.textContent = 'Pin plot';
                        pinBtn.setAttribute('aria-label', 'Pin ' + ((title.textContent || '').trim() || 'plot') + ' to scratchpad');

                        var downloadBtn = header.querySelector('.plot-download-btn, .shiny-download-link');
                        if (downloadBtn) {
                            header.insertBefore(pinBtn, downloadBtn);
                        } else {
                            header.appendChild(pinBtn);
                        }

                        card.__atlasPinReady = true;
                    });
                }

                function atlasTourState() {
                    if (!window.__atlasTourState) {
                        window.__atlasTourState = {
                            active: false,
                            stepIndex: 0,
                            markSeen: true
                        };
                    }
                    return window.__atlasTourState;
                }

                function atlasTourSteps() {
                    return [
                        {
                            title: 'Pick a source species and enter genes',
                            body: 'Open Camex or SATURN, choose the species that supplies your comparison genes, then stage one or more genes in the shared comparison panel.',
                            target: function() {
                                var tabs = Array.prototype.slice.call(document.querySelectorAll('.nav-tabs li a'));
                                return document.querySelector('.source-species-picker') ||
                                    document.getElementById('selected_genes') ||
                                    tabs.find(function(node) {
                                        return /Cross-species/.test(node.textContent || '');
                                    }) ||
                                    document.querySelector('.nav-tabs');
                            }
                        },
                        {
                            title: 'Apply the panel',
                            body: 'Click Apply comparison panel when the staged list looks right. That shared panel drives the cross-species comparisons, while species tabs keep their own native local panels.',
                            target: function() {
                                var tabs = Array.prototype.slice.call(document.querySelectorAll('.nav-tabs li a'));
                                return document.getElementById('apply_gene_selection') ||
                                    tabs.find(function(node) {
                                        return /Cross-species/.test(node.textContent || '');
                                    }) ||
                                    document.querySelector('.nav-tabs');
                            }
                        },
                        {
                            title: 'Browse a species tab',
                            body: 'The Medicago, Glycine, and Lotus tabs let you inspect native structure, markers, and local expression panels for each species.',
                            target: function() {
                                var tabs = Array.prototype.slice.call(document.querySelectorAll('.nav-tabs li a'));
                                return tabs.find(function(node) {
                                    var text = (node.textContent || '').trim();
                                    return text === 'Medicago truncatula' || text === 'Glycine max' || text === 'Lotus japonicus';
                                }) || document.querySelector('.nav-tabs');
                            }
                        },
                        {
                            title: 'Compare across species',
                            body: 'Use the Camex and SATURN tabs to inspect ortholog mappings and shared cross-species embeddings for the same panel.',
                            target: function() {
                                var tabs = Array.prototype.slice.call(document.querySelectorAll('.nav-tabs li a'));
                                return tabs.find(function(node) {
                                    return /Cross-species/.test(node.textContent || '');
                                }) || document.querySelector('.nav-tabs');
                            }
                        },
                        {
                            title: 'Download plots',
                            body: 'Use the plot download buttons in each card to save the current visualizations.',
                            target: function() {
                                return document.querySelector('.composite-dl-btn');
                            }
                        }
                    ];
                }

                function atlasClearTourTarget() {
                    var current = document.querySelector('.atlas-tour-target');
                    if (current) current.classList.remove('atlas-tour-target');
                }

                function atlasEnsureTourUi() {
                    var root = document.getElementById('atlas-tour-root');
                    if (root) return root;

                    root = document.createElement('div');
                    root.id = 'atlas-tour-root';
                    root.innerHTML = [
                        '<div class=\"atlas-tour-backdrop\"></div>',
                        '<div class=\"atlas-tour-card\" role=\"dialog\" aria-modal=\"true\" aria-live=\"polite\">',
                        '  <div class=\"atlas-tour-step\"></div>',
                        '  <h3 class=\"atlas-tour-title\"></h3>',
                        '  <p class=\"atlas-tour-body\"></p>',
                        '  <div class=\"atlas-tour-actions\">',
                        '    <button type=\"button\" class=\"atlas-tour-btn atlas-tour-btn-secondary\" data-tour-action=\"back\">Back</button>',
                        '    <button type=\"button\" class=\"atlas-tour-btn atlas-tour-btn-secondary\" data-tour-action=\"skip\">Skip</button>',
                        '    <button type=\"button\" class=\"atlas-tour-btn atlas-tour-btn-primary\" data-tour-action=\"next\">Next</button>',
                        '  </div>',
                        '</div>'
                    ].join('');
                    document.body.appendChild(root);

                    root.addEventListener('click', function(event) {
                        var action = event.target && event.target.getAttribute('data-tour-action');
                        if (!action) return;

                        var state = atlasTourState();
                        var lastStep = atlasTourSteps().length - 1;

                        if (action === 'back') {
                            state.stepIndex = Math.max(0, state.stepIndex - 1);
                            atlasRenderTour();
                            return;
                        }

                        if (action === 'skip') {
                            atlasStopTour(false);
                            return;
                        }

                        if (state.stepIndex >= lastStep) {
                            atlasStopTour(true);
                            return;
                        }

                        state.stepIndex += 1;
                        atlasRenderTour();
                    });

                    return root;
                }

                function atlasRenderTour() {
                    var state = atlasTourState();
                    var steps = atlasTourSteps();
                    var step = steps[state.stepIndex];
                    var root = atlasEnsureTourUi();
                    var card = root.querySelector('.atlas-tour-card');
                    var stepEl = root.querySelector('.atlas-tour-step');
                    var titleEl = root.querySelector('.atlas-tour-title');
                    var bodyEl = root.querySelector('.atlas-tour-body');
                    var backBtn = root.querySelector('[data-tour-action=\"back\"]');
                    var nextBtn = root.querySelector('[data-tour-action=\"next\"]');
                    var target = step.target ? step.target() : null;

                    atlasClearTourTarget();

                    if (target) {
                        target.classList.add('atlas-tour-target');
                        target.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    }

                    stepEl.textContent = 'Step ' + (state.stepIndex + 1) + ' of ' + steps.length;
                    titleEl.textContent = step.title;
                    bodyEl.textContent = step.body;
                    backBtn.disabled = state.stepIndex === 0;
                    nextBtn.textContent = state.stepIndex === steps.length - 1 ? 'Done' : 'Next';

                    root.className = 'atlas-tour-root is-active';

                    var defaultTop = Math.max(88, Math.round(window.innerHeight * 0.2));
                    var defaultLeft = Math.max(20, window.innerWidth - 380);

                    if (!target) {
                        card.style.top = defaultTop + 'px';
                        card.style.left = defaultLeft + 'px';
                        return;
                    }

                    var rect = target.getBoundingClientRect();
                    var cardWidth = 340;
                    var top = Math.max(88, Math.min(window.innerHeight - 220, rect.bottom + 18));
                    var left = rect.left;

                    if (left + cardWidth > window.innerWidth - 20) {
                        left = window.innerWidth - cardWidth - 20;
                    }
                    left = Math.max(20, left);

                    if (top > window.innerHeight - 220) {
                        top = Math.max(88, rect.top - 190);
                    }

                    card.style.top = Math.round(top) + 'px';
                    card.style.left = Math.round(left) + 'px';
                }

                function atlasStopTour(markSeen) {
                    var state = atlasTourState();
                    var root = document.getElementById('atlas-tour-root');

                    atlasClearTourTarget();

                    state.active = false;

                    if (root) {
                        root.className = 'atlas-tour-root';
                    }

                    if (markSeen && state.markSeen) {
                        try {
                            window.localStorage.setItem('atlas_tour_seen_v1', '1');
                        } catch (e) { /* localStorage blocked — ignore */ }
                    }
                }

                function atlasStartTour(markSeen) {
                    var overlay = document.getElementById('app-gate-overlay');
                    if (overlay && overlay.classList.contains('is-active')) {
                        window.__atlasPendingTourStart = !!markSeen;
                        return;
                    }

                    var state = atlasTourState();
                    state.active = true;
                    state.stepIndex = 0;
                    state.markSeen = markSeen !== false;

                    if (state.markSeen) {
                        try {
                            window.localStorage.setItem('atlas_tour_seen_v1', '1');
                        } catch (e) { /* localStorage blocked — ignore */ }
                    }

                    atlasRenderTour();
                }

                function initAtlasClient() {
                    if (!(window.Shiny && Shiny.setInputValue && Shiny.addCustomMessageHandler)) {
                        window.setTimeout(initAtlasClient, 250);
                        return;
                    }

                    if (window.__atlasClientInitialized) {
                        return;
                    }
                    window.__atlasClientInitialized = true;
                    window.__atlasPermalinkReportTimer = null;
                    window.__atlasPermalinkQuery = window.location.search || '';
                    window.__atlasRestoringPanelId = null;
                    window.__atlasBusyButtons = {};

                    atlasRenderPinboard();
                    atlasDecoratePlotCards();
                    [1500, 6000, 15000, 30000].forEach(function(delayMs) {
                        window.setTimeout(function() {
                            atlasDecoratePlotCards();
                            atlasReportVisiblePermalinkPanel();
                        }, delayMs);
                    });

                    try {
                        var seen = window.localStorage.getItem('atlas_tour_seen_v1');
                        if (!seen) {
                            window.setTimeout(function() {
                                atlasStartTour(true);
                            }, 900);
                        }
                    } catch (e) { /* localStorage blocked — skip tour */ }

                    var queryPayload = {};
                    var searchParams = new URLSearchParams(window.location.search);
                    searchParams.forEach(function(value, key) {
                        queryPayload[key] = value;
                    });
                    [1200, 5000, 15000, 30000].forEach(function(delayMs) {
                        window.setTimeout(function() {
                            var payload = Object.assign({}, queryPayload, { _nonce: Date.now() + delayMs });
                            Shiny.setInputValue('url_query_payload', payload, { priority: 'event' });
                            atlasReportVisiblePermalinkPanel();
                        }, delayMs);
                    });

                    var schedulePermalinkReport = function(delayMs) {
                        if (window.__atlasPermalinkReportTimer) {
                            window.clearTimeout(window.__atlasPermalinkReportTimer);
                        }
                        window.__atlasPermalinkReportTimer = window.setTimeout(function() {
                            atlasReportVisiblePermalinkPanel();
                        }, delayMs || 180);
                    };

                    var copyLink = document.getElementById('copy_permalink');
                    if (copyLink && !copyLink.__atlasBound) {
                        copyLink.__atlasBound = true;
                        copyLink.addEventListener('click', function(event) {
                            event.preventDefault();

                            var panelId = atlasVisiblePermalinkPanel();
                            var url = new URL(window.location.href);

                            if (window.__atlasPermalinkQuery && window.__atlasPermalinkQuery.length) {
                                url.search = window.__atlasPermalinkQuery;
                            }

                            if (panelId) {
                                url.searchParams.set('panel', panelId);
                            } else {
                                url.searchParams.delete('panel');
                            }

                            atlasReportVisiblePermalinkPanel(panelId);

                            navigator.clipboard.writeText(url.toString()).then(function() {
                                atlasShowToast('Copied permalink', 'success');
                            }).catch(function() {
                                atlasShowToast('Could not copy permalink', 'error');
                            });
                        });
                    }

                    var helpLink = document.getElementById('open_help_tour');
                    if (helpLink && !helpLink.__atlasTourBound) {
                        helpLink.__atlasTourBound = true;
                        helpLink.addEventListener('click', function(event) {
                            event.preventDefault();
                            atlasStartTour(false);
                        });
                    }

                    document.addEventListener('click', function(event) {
                        var tabLink = event.target && event.target.closest('.nav-tabs a');
                        if (tabLink) {
                            window.setTimeout(function() {
                                atlasReportVisiblePermalinkPanel();
                            }, 260);
                        }

                        var pinBtn = event.target && event.target.closest('.plot-card-pin-btn');
                        if (pinBtn) {
                            event.preventDefault();
                            atlasPinPlotCard(pinBtn.closest('.plot-card'));
                            return;
                        }

                        var removeBtn = event.target && event.target.closest('[data-pinboard-remove]');
                        if (removeBtn) {
                            event.preventDefault();
                            var state = atlasPinboardState();
                            var targetId = removeBtn.getAttribute('data-pinboard-remove');
                            state.items = state.items.filter(function(item) {
                                return item.id !== targetId;
                            });
                            atlasRenderPinboard();
                            return;
                        }

                        var clearBtn = event.target && event.target.closest('.atlas-pinboard-clear');
                        if (clearBtn) {
                            event.preventDefault();
                            atlasPinboardState().items = [];
                            atlasRenderPinboard();
                        }
                    });

                    window.addEventListener('scroll', function() {
                        schedulePermalinkReport(180);
                    }, { passive: true });
                    window.addEventListener('resize', function() {
                        schedulePermalinkReport(180);
                    });

                    document.addEventListener('shiny:idle', function() {
                        var busyButtons = window.__atlasBusyButtons || {};

                        Object.keys(busyButtons).forEach(function(buttonId) {
                            atlasSetButtonBusy(buttonId, false);
                        });

                        window.__atlasBusyButtons = {};
                        Shiny.setInputValue('atlas_client_idle', Date.now(), { priority: 'event' });
                    });

                    Shiny.addCustomMessageHandler('atlas_gate_unlock', function(msg) {
                        var overlay = document.getElementById('app-gate-overlay');
                        if (overlay) overlay.classList.remove('is-active');
                        if (window.__atlasPendingTourStart) {
                            window.__atlasPendingTourStart = false;
                            window.setTimeout(function() {
                                atlasStartTour(true);
                            }, 700);
                        }
                    });
                    Shiny.addCustomMessageHandler('atlas_restore_panel', function(msg) {
                        var panelId = msg && msg.panel ? msg.panel : null;
                        if (!panelId) return;

                        window.__atlasRestoringPanelId = panelId;
                        var attempts = 0;
                        var panelIsVisible = function(target) {
                            if (!target) return false;
                            var rect = target.getBoundingClientRect();
                            return rect.bottom > 80 && rect.top < window.innerHeight * 0.9;
                        };
                        var tryScroll = function() {
                            var target = atlasFindPermalinkPanel(panelId);
                            if (target && target.offsetParent !== null && panelIsVisible(target)) {
                                window.__atlasRestoringPanelId = null;
                                atlasReportVisiblePermalinkPanel(panelId);
                                return;
                            }

                            if (target && target.offsetParent !== null) {
                                target.scrollIntoView({ behavior: 'smooth', block: attempts < 2 ? 'start' : 'center' });
                                window.setTimeout(function() {
                                    var refreshedTarget = atlasFindPermalinkPanel(panelId);
                                    if (refreshedTarget && refreshedTarget.offsetParent !== null && panelIsVisible(refreshedTarget)) {
                                        window.__atlasRestoringPanelId = null;
                                        atlasReportVisiblePermalinkPanel(panelId);
                                        return;
                                    }
                                    if (attempts >= 40) {
                                        window.__atlasRestoringPanelId = null;
                                        atlasReportVisiblePermalinkPanel(panelId);
                                        return;
                                    }
                                    attempts += 1;
                                    window.setTimeout(tryScroll, 400);
                                }, 260);
                                return;
                            }

                            if (attempts >= 40) {
                                window.__atlasRestoringPanelId = null;
                                atlasReportVisiblePermalinkPanel(panelId);
                                return;
                            }
                            attempts += 1;
                            window.setTimeout(tryScroll, 400);
                        };

                        window.setTimeout(function() {
                            var initialTarget = atlasFindPermalinkPanel(panelId);
                            if (initialTarget && initialTarget.offsetParent !== null && panelIsVisible(initialTarget)) {
                                window.__atlasRestoringPanelId = null;
                                atlasReportVisiblePermalinkPanel(panelId);
                            } else {
                                tryScroll();
                            }
                        }, 250);
                    });
                    Shiny.addCustomMessageHandler('atlas_button_busy', function(msg) {
                        if (!msg || !msg.button_id) return;

                        var isBusy = !(msg.busy === false);
                        atlasSetButtonBusy(msg.button_id, isBusy, msg.label);

                        if (isBusy) {
                            window.__atlasBusyButtons[msg.button_id] = true;
                        } else if (window.__atlasBusyButtons) {
                            delete window.__atlasBusyButtons[msg.button_id];
                        }
                    });
                    Shiny.addCustomMessageHandler('atlas_set_permalink', function(msg) {
                        window.__atlasPermalinkQuery = (msg && msg.query) ? msg.query : '';
                    });
                }

                if (document.readyState === 'loading') {
                    document.addEventListener('DOMContentLoaded', initAtlasClient);
                } else {
                    initAtlasClient();
                }
            })();"
        ))
    ),
    uiOutput("app_gate_ui"),
    tags$div(id = "app-gate-overlay", class = if (atlas_access_required) "app-gate-overlay is-active" else "app-gate-overlay"),
        div(
            class = "app-shell",
            div(
                class = "busy-indicator",
                role = "status",
                `aria-live` = "polite",
                `aria-atomic` = "true",
                div(class = "busy-indicator-card",
                    div(class = "busy-indicator-spinner", `aria-hidden` = "true"),
                    div(class = "busy-indicator-copy",
                        div(class = "busy-indicator-text", "Updating...")
                    )
                )
            ),
        uiOutput("startup_data_check"),
        tags$header(
            class = "app-brand-shell",
            div(
                class = "app-brand-lockup",
                div(class = "app-brand-glyph", `aria-hidden` = "true", icon("leaf")),
                div(
                    class = "app-brand-text-wrap",
                    div(class = "app-brand-title", "Legume Root Nodule Symbiosis Atlas")
                )
            )
        ),
        div(
            class = "hero-copy hero-copy-wide",
            div(class = "hero-eyebrow", "Single-cell explorer"),
            h1("An scRNA-seq atlas for determinate and indeterminate nodules, with cross-species comparison"),
            p(
                class = "hero-text",
                HTML("Search by gene ID or common name to compare expression across <em>Medicago truncatula</em>, <em>Glycine max</em>, <em>Lotus japonicus</em>, and the cross-species integration.")
            ),
            div(
                class = "hero-badge-row",
                span(class = "hero-badge", "3 legume species"),
                span(class = "hero-badge", "within-species + cross-species views"),
                span(class = "hero-badge", "ortholog-aware expression plots"),
                span(class = "hero-badge hero-badge-warm", "pre-publication atlas")
            )
        ),
        uiOutput("atlas_summary_ui"),
        tags$main(
            id = "main-content",
            tabindex = "-1",
            div(
                class = "section-card",
                div(
                    class = "section-header",
                    div(class = "section-eyebrow", "Explorer settings"),
                    h2("Display and export"),
                    p("These settings apply across the atlas while you move between native species tabs and the cross-species integrations.")
                ),
                fluidRow(
                    column(
                        width = 4,
                        div(
                            class = "option-group",
                            selectInput(
                                inputId = "dl_format",
                                label = "Download format",
                                choices = c("SVG", "PNG", "PDF"),
                                selected = "PNG"
                            )
                        )
                    ),
                    column(
                        width = 4,
                        div(
                            class = "option-group",
                            selectInput(
                                inputId = "figure_preset",
                                label = "Figure style",
                                choices = c(
                                    "Exploratory" = "exploratory",
                                    "Presentation" = "presentation",
                                    "Publication" = "publication"
                                ),
                                selected = "publication"
                            )
                        )
                    ),
                    column(
                        width = 4,
                        div(
                            class = "option-group option-toggle",
                            checkboxInput(
                                inputId = "colorblind_safe",
                                label = "Colorblind-safe expression colors",
                                value = FALSE
                            )
                        )
                    )
                )
            ),
        do.call(
            tabsetPanel,
            c(
                list(
                    id = "main_tabs",
                    tabPanel(
                        title = "Overview",
                        value = "overview",
                        div(
                            class = "section-card overview-landing-card",
                            div(
                                class = "section-header",
                                div(class = "section-eyebrow", "Start here"),
                                h2("Choose how you want to explore the atlas"),
                                p("Use the species-specific tabs for native expression and marker discovery, then switch to a cross-species tab when you want orthology-aware comparison.")
                            ),
                            div(
                                class = "alert-stack overview-landing-stack",
                                notice_card(
                                    title = "Within-species tabs",
                                    body = "Build native local gene panels, inspect expression in each atlas, and use cluster markers without leaving the current species.",
                                    tone = "info"
                                ),
                                notice_card(
                                    title = "Cross-species tabs",
                                    body = "Build a shared comparison panel inside Camex or SATURN, choose the source species there, and review ortholog-aware mapping summaries directly in those tabs.",
                                    tone = "info"
                                )
                            )
                        ),
                        div(
                            class = "overview-feature-grid",
                            div(
                                class = "overview-feature-card",
                                div(class = "overview-feature-step", "1"),
                                h3("Start with a species"),
                                p("Use the Medicago, Glycine, or Lotus tabs when you already know the native gene IDs you want to inspect.")
                            ),
                            div(
                                class = "overview-feature-card",
                                div(class = "overview-feature-step", "2"),
                                h3("Stage a gene panel"),
                                p("Add genes by search, paste/import a list, or stage cluster markers. Plots refresh only after you click the Generate or Apply button.")
                            ),
                            div(
                                class = "overview-feature-card",
                                div(class = "overview-feature-step", "3"),
                                h3("Compare orthologs"),
                                p("Move to Camex or SATURN to see how selected source genes resolve across orthogroups and cross-species feature spaces.")
                            ),
                            div(
                                class = "overview-feature-card",
                                div(class = "overview-feature-step", "4"),
                                h3("Export with context"),
                                p("Download figures with atlas provenance, pin plots to the scratchpad, and copy a permalink for reproducible follow-up.")
                            )
                        ),
                        div(
                            class = "overview-release-note",
                            strong("Interpretation note: "),
                            "Orthogroup membership is a guide for exploration, not proof of one-to-one functional equivalence. One-to-many mappings and missing features are flagged in the cross-species tabs."
                        )
                    ),
                    species_tab_ui("medicago"),
                    species_tab_ui("glycine"),
                    species_tab_ui("lotus")
                ),
                lapply(cross_integration_keys, cross_tab_ui)
            )
        )
        ),
        div(
            id = "atlas-pinboard",
            class = "section-card atlas-pinboard",
            div(
                class = "section-header",
                div(class = "section-eyebrow", "Scratchpad"),
                h2("Pinned plots"),
                p("Pin up to four rendered plots from any tab to compare them side-by-side without losing the current view.")
            ),
            div(class = "atlas-pinboard-empty", "No plots pinned yet. Use “Pin plot” in any plot card to add the current view here."),
            div(class = "atlas-pinboard-actions",
                tags$button(type = "button", class = "plot-download-btn atlas-pinboard-clear", "Clear scratchpad")
            ),
            div(class = "atlas-pinboard-grid")
        ),
        tags$footer(
            class = "app-footer",
            div(
                class = "app-footer-row",
                div(
                    class = "app-footer-meta",
                    span(class = "app-footer-version", sprintf("Version %s", atlas_version)),
                    span(class = "app-footer-sep", aria_hidden = "true"),
                    span(class = "app-footer-updated", sprintf("Data updated %s", atlas_last_updated))
                ),
                div(
                    class = "app-footer-actions",
                    actionLink(inputId = "open_citation", label = "Cite this atlas", icon = icon("quote-right")),
                    actionLink(inputId = "copy_permalink", label = "Copy permalink", icon = icon("link")),
                    actionLink(inputId = "open_help_tour", label = "Help & walkthrough", icon = icon("circle-question"))
                )
            ),
            div(
                class = "app-footer-copy",
                "Single-cell atlas explorer for legume root nodule symbiosis.",
                strong(" 2026.")
            )
        )
    )
)

# ==============================================================================
# 3. SERVER LOGIC
# ==============================================================================

server <- function(input, output, session) {
    selection_notice <- reactiveVal(NULL)
    applied_selected_genes <- reactiveVal(character(0))
    comparison_source_species_state <- reactiveVal("medicago")
    comparison_staged_genes_state <- reactiveVal(character(0))
    gene_panel_busy_message <- reactiveVal(NULL)
    marker_job_messages <- reactiveVal(list())
    local_applied_gene_panels <- reactiveValues()
    local_gene_panel_busy_messages <- reactiveVal(setNames(as.list(rep(NA_character_, length(within_species_keys))), within_species_keys))
    local_selection_notices <- reactiveVal(setNames(as.list(rep(NA_character_, length(within_species_keys))), within_species_keys))
    gene_import_target <- reactiveVal("comparison")

    app_unlocked <- reactiveVal(!atlas_access_required)
    pending_url_genes <- reactiveVal(NULL)
    pending_url_local_species <- reactiveVal(NULL)
    pending_url_local_genes <- reactiveVal(NULL)
    pending_url_cluster <- reactiveVal(NULL)
    pending_url_panel <- reactiveVal(NULL)
    restoring_permalink_panel <- reactiveVal(NULL)
    permalink_panel_state <- reactiveVal(NULL)
    explore_mode <- reactive({
        input$explore_mode %||% "gene_driven"
    })

    output$app_gate_ui <- renderUI({
        if (isTRUE(app_unlocked())) return(NULL)
        div(
            class = "app-gate",
            role = "dialog",
            `aria-modal` = "true",
            `aria-labelledby` = "app-gate-title",
            div(
                class = "app-gate-card",
                div(class = "app-gate-icon", `aria-hidden` = "true", icon("leaf")),
                h2(id = "app-gate-title", class = "app-gate-title", "Pre-publication access"),
                p(
                    class = "app-gate-copy",
                    "This atlas is currently restricted while the accompanying manuscript is under peer review. Enter the access password shared with reviewers to continue."
                ),
                div(
                    class = "app-gate-form",
                    passwordInput("app_gate_password", label = "Password", placeholder = ""),
                    actionButton("app_gate_submit", "Unlock", class = "btn-primary"),
                    uiOutput("app_gate_error")
                )
            )
        )
    })

    observeEvent(input$app_gate_submit, {
        pw <- input$app_gate_password %||% ""
        if (identical(pw, atlas_access_password)) {
            app_unlocked(TRUE)
            session$sendCustomMessage("atlas_gate_unlock", TRUE)
        } else {
            output$app_gate_error <- renderUI({
                div(class = "app-gate-error", "Incorrect password.")
            })
        }
    })

    output$startup_data_check <- renderUI({
        missing <- check_missing_data_files()
        if (!length(missing)) {
            return(NULL)
        }

        bullets <- tags$ul(
            class = "startup-missing-list",
            lapply(missing, function(entry) {
                tags$li(tags$strong(entry$label), tags$code(entry$path))
            })
        )

        div(
            class = "alert-card warning startup-data-banner",
            div(class = "alert-title", "Atlas data is incomplete"),
            div(
                class = "alert-body",
                tagList(
                    tags$p("The following files were not found. Affected tabs will show an error until the data is in place:"),
                    bullets,
                    tags$p(
                        "If you are running a local copy, check that the data directory is mounted and that the ATLAS_DATA_DIR environment variable (if used) points to the right location."
                    )
                )
            )
        )
    })

    current_species_integration <- function(species_key) {
        input[[paste0(species_key, "_integration")]] %||% "ComBat_BBKNN"
    }

    local_panel_input_id <- function(species_key) {
        paste0(species_key, "_local_selected_genes")
    }

    local_panel_apply_button_id <- function(species_key) {
        paste0(species_key, "_apply_local_genes")
    }

    get_local_message <- function(store, species_key) {
        current <- store()
        value <- current[[species_key]]
        if (!length(value)) NA_character_ else value[[1]]
    }

    set_local_message <- function(store, species_key, message = NULL) {
        current <- isolate(store())
        current[[species_key]] <- if (is.null(message) || !nzchar(message)) NA_character_ else message
        store(current)
    }

    update_local_selected_genes_input <- function(species_key, choice_bundle, selected = character(0)) {
        updateSelectizeInput(
            session = session,
            inputId = local_panel_input_id(species_key),
            choices = choice_bundle$choices,
            selected = selected,
            server = TRUE
        )
    }

    local_staged_gene_values <- function(species_key) {
        genes <- input[[local_panel_input_id(species_key)]] %||% character(0)
        genes <- unique(as.character(genes))
        genes[!grepl("^__bulk__::", genes)]
    }

    local_applied_gene_values <- function(species_key) {
        genes <- local_applied_gene_panels[[species_key]] %||% character(0)
        unique(as.character(genes))
    }

    set_marker_job_message <- function(prefix, message = NULL) {
        current <- marker_job_messages()

        if (is.null(message) || !nzchar(message)) {
            current[[prefix]] <- NULL
        } else {
            current[[prefix]] <- message
        }

        marker_job_messages(current)
    }

    update_selected_genes_input <- function(choice_bundle, selected = character(0)) {
        comparison_staged_genes_state(unique(as.character(selected)))
        updateSelectizeInput(
            session = session,
            inputId = "selected_genes",
            choices = choice_bundle$choices,
            selected = selected,
            server = TRUE
        )
    }

    current_comparison_source_species <- reactive({
        comparison_source_species_state() %||% "medicago"
    })

    observeEvent(input$source_species, {
        selected_species <- input$source_species %||% "medicago"
        comparison_source_species_state(selected_species)
    }, ignoreInit = TRUE)

    observeEvent(input$selected_genes, {
        genes <- input$selected_genes %||% character(0)
        genes <- unique(as.character(genes))
        genes <- genes[!grepl("^__bulk__::", genes)]
        comparison_staged_genes_state(genes)
    }, ignoreInit = TRUE)

    staged_source_genes <- reactive({
        comparison_staged_genes_state()
    })

    selected_source_genes <- reactive({
        applied_selected_genes()
    })

    source_integration <- reactive({
        current_species_integration(current_comparison_source_species())
    })

    source_gene_catalog <- reactive({
        build_gene_catalog(
            current_comparison_source_species(),
            source_integration()
        )
    })

    observeEvent(
        source_gene_catalog(),
        {
            choice_bundle <- build_gene_choices(
                current_comparison_source_species(),
                source_integration()
            )

            current_selection <- isolate(staged_source_genes())
            current_applied_selection <- isolate(applied_selected_genes())
            valid_selection <- intersect(current_selection, choice_bundle$feature_ids)
            valid_applied_selection <- intersect(current_applied_selection, choice_bundle$feature_ids)
            dropped_genes <- unique(c(
                setdiff(current_selection, valid_selection),
                setdiff(current_applied_selection, valid_applied_selection)
            ))

            if (length(dropped_genes)) {
                selection_notice(
                    paste0(
                        "Dropped genes not present in the new source atlas: ",
                        compact_display_gene_list(current_comparison_source_species(), dropped_genes, limit = 6)
                    )
                )
                showNotification(
                    ui = tagList(
                        tags$strong("Some genes were removed from the panel."),
                        tags$br(),
                        sprintf(
                            "%d gene(s) are not present in the new source atlas: %s",
                            length(dropped_genes),
                            compact_display_gene_list(current_comparison_source_species(), dropped_genes, limit = 6)
                        )
                    ),
                    type = "warning",
                    duration = 10
                )
            } else {
                selection_notice(NULL)
            }

            applied_selected_genes(valid_applied_selection)

            update_selected_genes_input(
                choice_bundle = choice_bundle,
                selected = valid_selection
            )
        },
        ignoreInit = FALSE
    )

    observeEvent(
        input$apply_gene_selection,
        {
            staged_genes <- staged_source_genes()
            gene_panel_busy_message(
                if (length(staged_genes)) {
                    sprintf(
                        "Refreshing the shared comparison panel for %d staged gene(s)...",
                        length(staged_genes)
                    )
                } else {
                    "Clearing the shared comparison panel..."
                }
            )
            session$sendCustomMessage("atlas_button_busy", list(
                button_id = "apply_gene_selection",
                label = if (length(staged_genes)) "Applying..." else "Clearing..."
            ))
            applied_selected_genes(staged_genes)
            session$onFlushed(function() {
                gene_panel_busy_message(NULL)
                session$sendCustomMessage("atlas_button_busy", list(
                    button_id = "apply_gene_selection",
                    busy = FALSE
                ))
            }, once = TRUE)
        }
    )

    parse_pasted_gene_list <- function(raw) {
        if (is.null(raw) || !nzchar(raw)) return(character(0))
        tokens <- unlist(strsplit(raw, "[\\s,;\\t]+", perl = TRUE), use.names = FALSE)
        tokens <- trimws(tokens)
        tokens <- tokens[nzchar(tokens)]
        unique(tokens)
    }

    read_gene_csv <- function(path) {
        if (is.null(path) || !nzchar(path) || !file.exists(path)) return(character(0))
        sep <- if (grepl("\\.tsv$", path, ignore.case = TRUE)) "\t" else ","
        df <- tryCatch(
            utils::read.delim(path, sep = sep, stringsAsFactors = FALSE, check.names = FALSE, header = TRUE),
            error = function(e) NULL
        )
        if (is.null(df) || !ncol(df)) return(character(0))
        first_col <- df[[1]]
        first_col <- trimws(as.character(first_col))
        first_col <- first_col[nzchar(first_col)]
        unique(first_col)
    }

    observeEvent(input$open_citation, {
        showModal(modalDialog(
            title = "Cite this atlas",
            size = "m",
            easyClose = TRUE,
            tags$p("If you use data or visualizations from this explorer in your own work, please cite:"),
            tags$blockquote(class = "citation-block", atlas_citation_text),
            tags$p(
                class = "citation-meta",
                sprintf("App version %s — data snapshot %s.", atlas_version, atlas_last_updated)
            ),
            tags$p(
                class = "citation-meta",
                HTML("<em>A peer-reviewed publication describing this atlas is in preparation. This citation will be updated once the paper is published.</em>")
            ),
            footer = modalButton("Close")
        ))
    })

    url_restoring <- reactiveVal(TRUE)
    url_query_applied <- reactiveVal(FALSE)

    observeEvent(input$visible_permalink_panel, {
        payload <- input$visible_permalink_panel %||% list()
        panel_id <- payload[["panel"]] %||% ""
        panel_id <- trimws(as.character(panel_id[[1]] %||% panel_id))
        restoring_panel <- restoring_permalink_panel()

        if (!nzchar(panel_id)) {
            if (!is.null(restoring_panel) && nzchar(restoring_panel)) {
                return()
            }
            permalink_panel_state(NULL)
            return()
        }

        if (!is.null(restoring_panel) && nzchar(restoring_panel)) {
            if (!identical(panel_id, restoring_panel)) {
                return()
            }
            restoring_permalink_panel(NULL)
        }

        permalink_panel_state(panel_id)
    }, ignoreInit = TRUE)

    current_permalink_cluster <- reactive({
        current_tab <- input$main_tabs %||% "overview"

        if (current_tab %in% within_species_keys) {
            return(input[[paste0(current_tab, "_marker_cluster")]] %||% "")
        }

        if (current_tab %in% paste0("cross_", cross_integration_keys)) {
            return(input[[paste0(current_tab, "_marker_cluster")]] %||% "")
        }

        ""
    })

    build_query_string <- function(panel = permalink_panel_state()) {
        applied <- applied_selected_genes()
        species <- current_comparison_source_species()
        tab <- input$main_tabs %||% "overview"
        cb <- if (isTRUE(input$colorblind_safe)) "1" else "0"
        mode <- explore_mode()

        parts <- c(
            sprintf("mode=%s", utils::URLencode(mode, reserved = TRUE)),
            sprintf("species=%s", utils::URLencode(species, reserved = TRUE)),
            sprintf("tab=%s", utils::URLencode(tab, reserved = TRUE)),
            sprintf("cb=%s", cb)
        )

        for (species_key in within_species_keys) {
            integration_method <- current_species_integration(species_key)

            if (!is.null(integration_method) && nzchar(integration_method)) {
                parts <- c(
                    parts,
                    sprintf(
                        "%s_integration=%s",
                        species_key,
                        utils::URLencode(integration_method, reserved = TRUE)
                    )
                )
            }
        }

        cluster_id <- current_permalink_cluster()
        if (!is.null(cluster_id) && nzchar(cluster_id)) {
            parts <- c(parts, sprintf("cluster=%s", utils::URLencode(cluster_id, reserved = TRUE)))
        }

        if (!is.null(panel) && nzchar(panel)) {
            parts <- c(parts, sprintf("panel=%s", utils::URLencode(panel, reserved = TRUE)))
        }

        if (length(applied)) {
            parts <- c(parts, sprintf("genes=%s", utils::URLencode(paste(applied, collapse = ","), reserved = TRUE)))
        }

        if (tab %in% within_species_keys) {
            local_genes <- local_applied_gene_values(tab)

            if (length(local_genes)) {
                parts <- c(
                    parts,
                    sprintf("local_species=%s", utils::URLencode(tab, reserved = TRUE)),
                    sprintf("local_genes=%s", utils::URLencode(paste(local_genes, collapse = ","), reserved = TRUE))
                )
            }
        }

        paste0("?", paste(parts, collapse = "&"))
    }

    observeEvent(input$url_query_payload, {
        if (isTRUE(url_query_applied())) {
            return()
        }

        query <- input$url_query_payload %||% list()

        if (!length(query)) {
            url_query_applied(TRUE)
            url_restoring(FALSE)
            return()
        }

        mode <- query[["mode"]]
        if (!is.null(mode) && nzchar(mode) && mode %in% c("gene_driven", "cluster_driven")) {
            updateRadioButtons(session, "explore_mode", selected = mode)
        }

        species <- query[["species"]]
        if (!is.null(species) && nzchar(species) && (species %in% names(species_choices) || species %in% unname(species_choices))) {
            comparison_source_species_state(species)
            shinyWidgets::updatePickerInput(session, "source_species", selected = species)
        }

        for (species_key in within_species_keys) {
            integration_param <- query[[paste0(species_key, "_integration")]]
            if (!is.null(integration_param) && nzchar(integration_param) && integration_param %in% unname(integration_choices)) {
                updateSelectInput(
                    session,
                    inputId = paste0(species_key, "_integration"),
                    selected = integration_param
                )
            }
        }

        tab <- query[["tab"]]
        if (!is.null(tab) && nzchar(tab)) {
            updateTabsetPanel(session, "main_tabs", selected = tab)
        }

        genes_param <- query[["genes"]]
        if (!is.null(genes_param) && nzchar(genes_param)) {
            pending_url_genes(strsplit(genes_param, ",", fixed = TRUE)[[1]])
        }

        local_species <- query[["local_species"]]
        local_genes <- query[["local_genes"]]
        if (!is.null(local_species) && nzchar(local_species) && local_species %in% within_species_keys &&
            !is.null(local_genes) && nzchar(local_genes)) {
            pending_url_local_species(local_species)
            pending_url_local_genes(strsplit(local_genes, ",", fixed = TRUE)[[1]])
        }

        cb <- query[["cb"]]
        if (!is.null(cb) && nzchar(cb)) {
            updateCheckboxInput(session, "colorblind_safe", value = identical(cb, "1"))
        }

        cluster_param <- query[["cluster"]]
        if (!is.null(cluster_param) && nzchar(cluster_param)) {
            pending_url_cluster(cluster_param)
        }

        panel_param <- query[["panel"]]
        if (!is.null(panel_param) && nzchar(panel_param)) {
            pending_url_panel(panel_param)
            restoring_permalink_panel(panel_param)
            permalink_panel_state(panel_param)
        }

        url_query_applied(TRUE)
        url_restoring(FALSE)
    }, ignoreInit = FALSE)

    observe({
        source_gene_catalog()
        pending <- pending_url_genes()
        if (is.null(pending) || !length(pending)) return()

        bundle <- build_gene_choices(
            current_comparison_source_species(),
            source_integration()
        )
        valid_ids <- bundle$feature_ids %||% character(0)
        pending_valid <- intersect(pending, valid_ids)
        if (length(pending_valid)) {
            update_selected_genes_input(
                choice_bundle = bundle,
                selected = pending_valid
            )
            applied_selected_genes(pending_valid)
        }
        pending_url_genes(NULL)
    })

    observeEvent(input$selected_genes_bulk_query, {
        payload <- input$selected_genes_bulk_query %||% list()
        query <- trimws(as.character(payload$query %||% ""))

        if (!nzchar(query)) {
            return()
        }

        source_species <- current_comparison_source_species()
        bundle <- build_gene_choices(source_species, source_integration())
        matched_ids <- match_gene_choices(bundle, query)
        current_selection <- isolate(staged_source_genes())
        new_matches <- setdiff(matched_ids, current_selection)
        updated_selection <- unique(c(current_selection, matched_ids))

        update_selected_genes_input(
            choice_bundle = bundle,
            selected = updated_selection
        )

        if (!length(matched_ids)) {
            showNotification(
                sprintf("No genes in %s match \"%s\".", species_label(source_species), query),
                type = "warning",
                duration = 6
            )
            return()
        }

        if (!length(new_matches)) {
            showNotification(
                sprintf(
                    "All %d gene(s) matching \"%s\" are already staged.",
                    length(matched_ids),
                    query
                ),
                type = "message",
                duration = 6
            )
            return()
        }

        showNotification(
            sprintf(
                "Added %d gene(s) matching \"%s\" to the staged panel.",
                length(new_matches),
                query
            ),
            type = "message",
            duration = 6
        )
    }, ignoreInit = TRUE)

    observeEvent(input$atlas_client_idle, {
        gene_panel_busy_message(NULL)
        marker_job_messages(list())
        local_gene_panel_busy_messages(
            setNames(as.list(rep(NA_character_, length(within_species_keys))), within_species_keys)
        )
    }, ignoreInit = TRUE)

    observe({
        if (isTRUE(url_restoring())) return()
        session$sendCustomMessage("atlas_set_permalink", list(
            query = build_query_string()
        ))
    })

    observe({
        panel_id <- pending_url_panel()

        if (isTRUE(url_restoring()) || is.null(panel_id) || !nzchar(panel_id) || !isTRUE(app_unlocked())) {
            return()
        }

        session$sendCustomMessage("atlas_restore_panel", list(panel = panel_id))
        pending_url_panel(NULL)
    })

    gene_import_context <- function(target_key = "comparison") {
        if (identical(target_key, "comparison")) {
            list(
                target_key = "comparison",
                species_key = current_comparison_source_species(),
                integration_method = source_integration(),
                panel_label = "cross-species comparison panel",
                species_label = species_label(current_comparison_source_species()),
                replace_label = "Replace current comparison panel instead of appending"
            )
        } else {
            list(
                target_key = target_key,
                species_key = target_key,
                integration_method = current_species_integration(target_key),
                panel_label = paste(species_label(target_key), "local panel"),
                species_label = species_label(target_key),
                replace_label = paste("Replace current", species_label(target_key), "local panel instead of appending")
            )
        }
    }

    show_gene_import_modal <- function(target_key = "comparison") {
        context <- gene_import_context(target_key)
        gene_import_target(target_key)

        showModal(modalDialog(
            title = "Import gene list",
            size = "m",
            easyClose = TRUE,
            tags$p(
                class = "modal-hint",
                sprintf(
                    "Paste or upload a list of gene IDs (one per line, or separated by commas / tabs / semicolons). Genes are matched against the current %s catalog for the %s; unrecognized genes are dropped.",
                    context$species_label,
                    context$panel_label
                )
            ),
            tags$label(`for` = "gene_import_paste", "Paste gene IDs"),
            tags$textarea(
                id = "gene_import_paste",
                class = "form-control gene-import-textarea",
                rows = 8,
                placeholder = "e.g.\nMedtr7g029290\nMedtr4g104370\n..."
            ),
            tags$div(class = "modal-divider", "or"),
            fileInput(
                inputId = "gene_import_file",
                label = "Upload CSV / TSV (first column is the gene ID)",
                accept = c(".csv", ".tsv", ".txt"),
                buttonLabel = "Choose file",
                placeholder = "No file selected"
            ),
            checkboxInput(
                inputId = "gene_import_replace",
                label = context$replace_label,
                value = FALSE
            ),
            footer = tagList(
                modalButton("Cancel"),
                actionButton("confirm_gene_import", "Import", class = "btn-primary")
            )
        ))
    }

    observeEvent(input$open_gene_import, {
        show_gene_import_modal("comparison")
    })

    observeEvent(input$confirm_gene_import, {
        pasted <- parse_pasted_gene_list(input$gene_import_paste)
        file_df <- input$gene_import_file
        from_file <- if (!is.null(file_df) && nrow(file_df)) read_gene_csv(file_df$datapath[1]) else character(0)
        imported <- unique(c(pasted, from_file))

        if (!length(imported)) {
            showNotification("No gene IDs detected in the pasted text or uploaded file.", type = "warning", duration = 6)
            return()
        }

        target_key <- gene_import_target() %||% "comparison"
        context <- gene_import_context(target_key)
        bundle <- build_gene_choices(
            context$species_key,
            context$integration_method
        )
        valid_ids <- bundle$feature_ids

        if (is.null(valid_ids)) valid_ids <- character(0)

        matches <- imported[imported %in% valid_ids]

        lower_map <- setNames(valid_ids, tolower(valid_ids))
        case_insensitive_hits <- lower_map[tolower(imported[!(imported %in% valid_ids)])]
        case_insensitive_hits <- unname(case_insensitive_hits[!is.na(case_insensitive_hits)])
        matches <- unique(c(matches, case_insensitive_hits))

        unmatched <- setdiff(imported, matches)
        unmatched <- unmatched[!(tolower(unmatched) %in% tolower(matches))]

        current <- if (identical(target_key, "comparison")) {
            isolate(staged_source_genes())
        } else {
            isolate(local_staged_gene_values(target_key))
        }
        new_selection <- if (isTRUE(input$gene_import_replace)) matches else unique(c(current, matches))

        if (identical(target_key, "comparison")) {
            update_selected_genes_input(
                choice_bundle = bundle,
                selected = new_selection
            )
        } else {
            update_local_selected_genes_input(
                species_key = target_key,
                choice_bundle = bundle,
                selected = new_selection
            )
        }

        removeModal()

        summary_msg <- sprintf(
            "Imported %d gene(s) into the %s%s. %d not recognized in the current %s catalog%s.",
            length(matches),
            context$panel_label,
            if (isTRUE(input$gene_import_replace)) " (replaced current selection)" else " (appended)",
            length(unmatched),
            context$species_label,
            if (length(unmatched)) paste0(": ", compact_display_gene_list(context$species_key, unmatched, limit = 6)) else ""
        )
        showNotification(summary_msg, type = if (length(unmatched)) "warning" else "message", duration = 8)
    })

    source_orthogroup_status <- reactive({
        genes <- selected_source_genes()

        if (!length(genes)) {
            return(list(
                with_orthogroup = character(0),
                without_orthogroup = character(0)
            ))
        }

        source_orthogroups <- resolve_source_orthogroups(
            current_comparison_source_species(),
            genes
        )

        genes_with_orthogroup <- source_orthogroups %>%
            filter(!is.na(orthogroup)) %>%
            distinct(source_gene) %>%
            pull(source_gene)

        list(
            with_orthogroup = genes_with_orthogroup,
            without_orthogroup = setdiff(genes, genes_with_orthogroup)
        )
    })

    overview_cross_resolutions <- reactive({
        source_species <- current_comparison_source_species()
        genes <- selected_source_genes()
        setNames(
            lapply(cross_integration_keys, function(cross_key) {
                resolve_cross_integration_mapping(
                    source_species = source_species,
                    source_genes = genes,
                    cross_key = cross_key
                )
            }),
            cross_integration_keys
        )
    }) %>% bindCache(
        current_comparison_source_species(),
        selected_source_genes(),
        cache = "app"
    )

    output$gene_selection_status <- renderText({
        staged_genes <- staged_source_genes()
        applied_genes <- selected_source_genes()
        n_genes <- length(applied_genes)
        note <- selection_notice()
        busy_message <- gene_panel_busy_message()

        if (!is.null(busy_message) && nzchar(busy_message)) {
            return(busy_message)
        }

        if (!same_gene_selection(staged_genes, applied_genes)) {
            staged_message <- if (length(staged_genes)) {
                paste0(
                    length(staged_genes),
                    " gene(s) staged. Click Apply comparison panel to refresh the cross-species views."
                )
            } else {
                "No genes staged. Click Apply comparison panel to clear the cross-species views."
            }

            if (!is.null(note)) {
                return(paste(staged_message, note))
            }

            return(staged_message)
        }

        if (n_genes == 0) {
            if (!is.null(note)) {
                return(paste("No comparison genes selected yet.", note))
            }

            return("No comparison genes selected yet.")
        }

        base_message <- paste0(
            n_genes,
            " gene(s) in the shared cross-species comparison panel."
        )

        if (!is.null(note)) {
            paste(base_message, note)
        } else {
            base_message
        }
    })

    output$atlas_summary_ui <- renderUI({
        within_cards <- lapply(within_species_keys, function(species_key) {
            summary <- get_within_dataset_summary(
                species_key,
                current_species_integration(species_key)
            )

            summary_strip_tile(species_label_tag(species_key), summary)
        })

        cross_summary <- get_cross_dataset_summary(cross_integration_keys[[1]])

        cross_card <- summary_strip_tile(
            tags$span(style = "text-transform: none; letter-spacing: normal;", "Cross-species"),
            cross_summary
        )

        div(
            class = "atlas-summary-strip",
            tagList(within_cards, cross_card)
        )
    })

    overview_mapping_table <- reactive({
        genes <- selected_source_genes()

        if (!length(genes)) {
            return(tibble())
        }

        source_species <- current_comparison_source_species()
        source_orthogroups <- resolve_source_orthogroups(source_species, genes) %>%
            distinct(source_gene, orthogroup, .keep_all = TRUE) %>%
            arrange(match(source_gene, genes))

        source_orthogroups %>%
            mutate(
                `Source species` = species_label(source_species),
                `Source gene` = display_gene_labels(source_species, source_gene),
                Orthogroup = ifelse(is.na(orthogroup), "No orthogroup", orthogroup),
                Medicago = map_chr(orthogroup, ~ compact_display_gene_list("medicago", get_orthogroup_members(.x, "medicago"))),
                `Glycine max` = map_chr(orthogroup, ~ compact_display_gene_list("glycine", get_orthogroup_members(.x, "glycine"))),
                `Lotus japonicus` = map_chr(orthogroup, ~ compact_display_gene_list("lotus", get_orthogroup_members(.x, "lotus")))
            ) %>%
            select(
                `Source species`,
                `Source gene`,
                Orthogroup,
                Medicago,
                `Glycine max`,
                `Lotus japonicus`
            )
    })

    output$overview_mapping_table_ui <- renderUI({
        mapping_tbl <- overview_mapping_table()

        if (!nrow(mapping_tbl)) {
            return(div(
                class = "summary-placeholder",
                "Select one or more source-species genes to inspect the ortholog mapping summary."
            ))
        }

        html_summary_table(mapping_tbl)
    })

    output$dl_overview_mapping_table <- downloadHandler(
        filename = function() paste0("ortholog_mapping_", Sys.Date(), ".csv"),
        content = function(file) {
            export_tbl <- add_export_provenance_columns(
                overview_mapping_table(),
                tab_label = "Overview",
                integration_label = "Cross-atlas overview",
                extra = list(atlas_export_type = "ortholog_mapping_table")
            )

            write.csv(export_tbl, file, row.names = FALSE)
        }
    )

    output$overview_alerts_ui <- renderUI({
        genes <- selected_source_genes()

        if (!length(genes)) {
            return(
                div(
                    class = "alert-stack",
                    notice_card(
                        title = "Ready for selection",
                        body = "Build a shared comparison panel here, or promote a local species panel from a within-species tab when you are ready for cross-species comparison.",
                        tone = "info"
                    )
                )
            )
        }

        source_species <- current_comparison_source_species()
        orthogroup_status <- source_orthogroup_status()
        other_species <- setdiff(within_species_keys, source_species)

        cards <- list(
            notice_card(
                title = "How to read ortholog mappings",
                body = "Orthogroup membership is not proof of one-to-one equivalence. One-to-many mappings can reflect paralogs, and missing members or features can come from orthology-table or atlas-coverage limits rather than true biological absence.",
                tone = "info"
            )
        )

        if (length(orthogroup_status$without_orthogroup)) {
            cards <- append(cards, list(
                notice_card(
                    title = "Selected genes without orthogroups",
                    body = paste(
                        compact_display_gene_list(source_species, orthogroup_status$without_orthogroup, limit = 8),
                        "These genes are absent from the current orthogroup table, so cross-species comparisons cannot be interpreted for them here."
                    ),
                    tone = "warning"
                )
            ))
        }

        for (target_species in other_species) {
            resolution <- resolve_target_mapping(
                source_species = source_species,
                source_genes = genes,
                target_species = target_species,
                integration_method = current_species_integration(target_species),
                cross_space = FALSE
            )

            if (length(resolution$no_target_members)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(species_label(target_species), "orthogroup gaps"),
                        body = paste(
                            "Orthogroups were found, but they do not contain target-species members for:",
                            compact_display_gene_list(source_species, resolution$no_target_members, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }

            if (length(resolution$missing_features)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(species_label(target_species), "feature gaps"),
                        body = paste(
                            "Orthologs were found, but no matching features are present in the selected atlas for:",
                            compact_display_gene_list(source_species, resolution$missing_features, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }

        }

        cross_resolutions <- overview_cross_resolutions()

        for (cross_key in cross_integration_keys) {
            cross_res <- cross_resolutions[[cross_key]]
            cross_label <- cross_integration_label(cross_key)

            if (length(cross_res$no_target_members)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(cross_label, "orthogroup gaps"),
                        body = paste(
                            "The orthogroups below lack usable target members in the",
                            cross_label,
                            "space for:",
                            compact_display_gene_list(source_species, cross_res$no_target_members, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }

            if (length(cross_res$missing_features)) {
                cards <- append(cards, list(
                    notice_card(
                        title = paste(cross_label, "feature gaps"),
                        body = paste(
                            "Orthologs were resolved but no matching features exist in the",
                            cross_label,
                            "object for:",
                            compact_display_gene_list(source_species, cross_res$missing_features, limit = 8)
                        ),
                        tone = "warning"
                    )
                ))
            }
        }

        if (!length(cards)) {
            cards <- list(
                notice_card(
                    title = "Mappings look consistent",
                    body = "The current gene panel resolves cleanly across the species-specific atlases and both cross-species integrations.",
                    tone = "info"
                )
            )
        }

        div(class = "alert-stack", tagList(cards))
    })

    get_ext <- function() {
        tolower(input$dl_format %||% "png")
    }

    figure_preset <- reactive({
        input$figure_preset %||% "publication"
    })

    observeEvent(figure_preset(), {
        recommended_format <- switch(
            figure_preset(),
            exploratory = "PNG",
            presentation = "SVG",
            publication = "PDF",
            "PDF"
        )

        updateSelectInput(session, "dl_format", selected = recommended_format)
    }, ignoreInit = TRUE)

    selected_gene_export_summary <- function(limit = 4L) {
        source_species <- input$source_species %||% "medicago"
        gene_labels <- display_gene_labels(
            source_species,
            selected_source_genes(),
            include_gene_id_with_common = FALSE
        )

        if (!length(gene_labels)) {
            return("none")
        }

        compact_gene_list(gene_labels, limit = limit)
    }

    export_provenance_fields <- function(tab_label, integration_label = NULL, extra = list()) {
        base_fields <- list(
            atlas_version = atlas_version,
            atlas_exported_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
            atlas_source_species = species_label(current_comparison_source_species()),
            atlas_tab = tab_label %||% (input$main_tabs %||% "overview"),
            atlas_integration = integration_label %||% "",
            atlas_selected_genes = selected_gene_export_summary()
        )

        c(base_fields, extra)
    }

    add_export_provenance_columns <- function(tbl, tab_label, integration_label = NULL, extra = list()) {
        tbl <- as_tibble(tbl)
        provenance_tbl <- as_tibble(export_provenance_fields(
            tab_label = tab_label,
            integration_label = integration_label,
            extra = extra
        ))

        if (!nrow(tbl)) {
            return(bind_cols(provenance_tbl[rep(1, 0), , drop = FALSE], tbl))
        }

        bind_cols(provenance_tbl[rep(1, nrow(tbl)), , drop = FALSE], tbl)
    }

    export_plot_caption <- function(tab_label, integration_label = NULL, extra = list()) {
        fields <- export_provenance_fields(tab_label, integration_label, extra = extra)
        extra_fields <- extra

        caption_parts <- c(
            paste0("Atlas ", fields$atlas_version),
            paste0("Exported ", fields$atlas_exported_at),
            paste0("Tab: ", fields$atlas_tab),
            paste0("Source species: ", fields$atlas_source_species),
            if (nzchar(fields$atlas_integration)) paste0("Integration: ", fields$atlas_integration) else NULL,
            paste0("Genes: ", fields$atlas_selected_genes),
            if (length(extra_fields)) {
                paste(
                    vapply(names(extra_fields), function(name) {
                        paste0(gsub("^atlas_", "", name), ": ", extra_fields[[name]])
                    }, character(1)),
                    collapse = " | "
                )
            } else {
                NULL
            }
        )

        paste(caption_parts[nzchar(caption_parts)], collapse = " | ")
    }

    add_plot_export_provenance <- function(plot_obj, tab_label, integration_label = NULL, extra = list()) {
        caption_text <- export_plot_caption(
            tab_label = tab_label,
            integration_label = integration_label,
            extra = extra
        )

        if (!nzchar(caption_text)) {
            return(plot_obj)
        }

        tryCatch(
            plot_obj + labs(caption = caption_text) + theme(
                plot.caption = element_text(
                    size = 8.5,
                    colour = app_palette["muted"],
                    hjust = 0
                )
            ),
            error = function(e) plot_obj
        )
    }

    apply_export_figure_preset <- function(plot_obj) {
        preset <- figure_preset()
        preset_theme <- app_plot_theme(preset = preset)

        tryCatch(
            plot_obj & preset_theme,
            error = function(e) plot_obj + preset_theme
        )
    }

    save_ggplot <- function(file, plot_obj, width, height, tab_label = NULL, integration_label = NULL, extra = list()) {
        preset_cfg <- figure_preset_config(figure_preset())
        plot_obj <- apply_export_figure_preset(plot_obj)
        plot_obj <- add_plot_export_provenance(
            plot_obj,
            tab_label = tab_label,
            integration_label = integration_label,
            extra = extra
        )

        ggsave(
            filename = file,
            plot = plot_obj,
            device = get_ext(),
            width = width * preset_cfg$width_scale,
            height = height * preset_cfg$height_scale,
            dpi = 300,
            limitsize = FALSE
        )
    }

    register_species_tab <- function(species_key) {
        local({
            prefix <- species_key
            tab_label <- species_label(species_key)
            clustering_columns <- unique(c(
                paste0(setdiff(distribution_cluster_columns, "cluster_label"), "_label"),
                "cluster_label",
                distribution_cluster_columns
            ))
            tab_dataset_key <- reactive({
                paste(species_key, current_species_integration(species_key), sep = "_")
            })
            tab_is_active <- reactive({
                identical(input$main_tabs %||% "overview", species_key)
            })
            tab_local_choice_bundle <- reactive({
                build_gene_choices(species_key, current_species_integration(species_key))
            })
            tab_local_staged_genes <- reactive({
                local_staged_gene_values(species_key)
            })
            tab_local_applied_genes <- reactive({
                local_applied_gene_values(species_key)
            })
            empty_local_selection_warning <- sprintf(
                "Select at least one gene to generate the expression plots for %s.",
                tab_label
            )

            observeEvent(tab_local_choice_bundle(), {
                choice_bundle <- tab_local_choice_bundle()
                current_selection <- isolate(tab_local_staged_genes())
                current_applied_selection <- isolate(tab_local_applied_genes())
                valid_selection <- intersect(current_selection, choice_bundle$feature_ids %||% character(0))
                valid_applied_selection <- intersect(current_applied_selection, choice_bundle$feature_ids %||% character(0))
                dropped_genes <- unique(c(
                    setdiff(current_selection, valid_selection),
                    setdiff(current_applied_selection, valid_applied_selection)
                ))

                if (length(dropped_genes)) {
                    set_local_message(
                        local_selection_notices,
                        species_key,
                        paste0(
                            "Dropped genes not present in the current ",
                            tab_label,
                            " atlas: ",
                            compact_display_gene_list(species_key, dropped_genes, limit = 6)
                        )
                    )
                    showNotification(
                        ui = tagList(
                            tags$strong(paste("Some", tab_label, "genes were removed from the local panel.")),
                            tags$br(),
                            compact_display_gene_list(species_key, dropped_genes, limit = 6)
                        ),
                        type = "warning",
                        duration = 8
                    )
                } else {
                    set_local_message(local_selection_notices, species_key, NULL)
                }

                local_applied_gene_panels[[species_key]] <- valid_applied_selection
                update_local_selected_genes_input(
                    species_key = species_key,
                    choice_bundle = choice_bundle,
                    selected = valid_selection
                )

                pending_species <- pending_url_local_species()
                pending_genes <- pending_url_local_genes()

                if (identical(pending_species, species_key) && length(pending_genes)) {
                    pending_valid <- intersect(pending_genes, choice_bundle$feature_ids %||% character(0))
                    update_local_selected_genes_input(
                        species_key = species_key,
                        choice_bundle = choice_bundle,
                        selected = pending_valid
                    )
                    local_applied_gene_panels[[species_key]] <- pending_valid
                    pending_url_local_species(NULL)
                    pending_url_local_genes(NULL)
                }
            }, ignoreInit = FALSE)

            observeEvent(tab_local_staged_genes(), {
                current_note <- get_local_message(local_selection_notices, species_key)
                if (length(tab_local_staged_genes()) > 0 && identical(current_note, empty_local_selection_warning)) {
                    set_local_message(local_selection_notices, species_key, NULL)
                }
            }, ignoreInit = TRUE)

            output[[paste0(prefix, "_local_gene_selection_status")]] <- renderText({
                staged_genes <- tab_local_staged_genes()
                applied_genes <- tab_local_applied_genes()
                busy_message <- get_local_message(local_gene_panel_busy_messages, species_key)
                note <- get_local_message(local_selection_notices, species_key)

                if (!is.na(busy_message) && nzchar(busy_message)) {
                    return(busy_message)
                }

                if (!same_gene_selection(staged_genes, applied_genes)) {
                    staged_message <- if (length(staged_genes)) {
                        paste0(length(staged_genes), " ", tab_label, " gene(s) staged. Click Generate the expression plots to refresh this tab.")
                    } else {
                        paste0("No ", tab_label, " genes staged. Select at least one gene, then click Generate the expression plots.")
                    }

                    if (!is.na(note) && nzchar(note)) {
                        return(paste(staged_message, note))
                    }

                    return(staged_message)
                }

                if (!length(applied_genes)) {
                    if (!is.na(note) && nzchar(note)) {
                        return(paste("No local genes selected yet.", note))
                    }

                    return("No local genes selected yet.")
                }

                base_message <- paste0(length(applied_genes), " ", tab_label, " gene(s) in the local panel.")

                if (!is.na(note) && nzchar(note)) {
                    paste(base_message, note)
                } else {
                    base_message
                }
            })

            observeEvent(input[[paste0(prefix, "_apply_local_genes")]], {
                staged_genes <- tab_local_staged_genes()

                if (!length(staged_genes)) {
                    set_local_message(local_gene_panel_busy_messages, species_key, NULL)
                    set_local_message(local_selection_notices, species_key, empty_local_selection_warning)
                    showNotification(
                        empty_local_selection_warning,
                        type = "warning",
                        duration = 6
                    )
                    return()
                }

                set_local_message(
                    local_gene_panel_busy_messages,
                    species_key,
                    sprintf("Refreshing %s plots for %d local gene(s)...", tab_label, length(staged_genes))
                )

                session$sendCustomMessage("atlas_button_busy", list(
                    button_id = local_panel_apply_button_id(species_key),
                    label = "Applying..."
                ))
                local_applied_gene_panels[[species_key]] <- staged_genes
                session$onFlushed(function() {
                    set_local_message(local_gene_panel_busy_messages, species_key, NULL)
                    session$sendCustomMessage("atlas_button_busy", list(
                        button_id = local_panel_apply_button_id(species_key),
                        busy = FALSE
                    ))
                }, once = TRUE)
            }, ignoreInit = TRUE)

            observeEvent(input[[paste0(prefix, "_open_gene_import")]], {
                show_gene_import_modal(species_key)
            }, ignoreInit = TRUE)

            tab_object <- reactive({
                req(tab_is_active())
                get_within_object(
                    species_key,
                    current_species_integration(species_key)
                )
            })

            tab_resolution <- reactive({
                req(tab_is_active())
                enforce_plot_feature_limit(resolve_target_mapping(
                    source_species = species_key,
                    source_genes = tab_local_applied_genes(),
                    target_species = species_key,
                    integration_method = current_species_integration(species_key),
                    cross_space = FALSE
                ))
            }) %>% bindCache(
                species_key,
                tab_local_applied_genes(),
                current_species_integration(species_key),
                atlas_plot_feature_limit,
                cache = "app"
            )

            expression_prompt <- sprintf(
                "Apply a local gene panel first to populate the expression panels for %s.",
                tab_label
            )
            validate_local_expression_ready <- function() {
                validate(need(length(tab_local_applied_genes()) > 0, expression_prompt))
            }

            expression_resolution <- reactive({
                validate(
                    need(
                        length(tab_local_applied_genes()) > 0,
                        expression_prompt
                    )
                )

                tab_resolution()
            })

            tab_group_choices <- reactive({
                get_cached_ui_choices(tab_dataset_key(), "within_group") %||%
                    within_group_choices(tab_object())
            })
            tab_distribution_split_choices <- reactive({
                get_cached_ui_choices(tab_dataset_key(), "within_distribution_split") %||%
                    within_distribution_split_choices(tab_object())
            })
            tab_feature_split_choices <- reactive({
                get_cached_ui_choices(tab_dataset_key(), "within_feature_split") %||%
                    within_feature_split_choices(tab_object())
            })
            tab_composition_choices <- reactive({
                get_cached_ui_choices(tab_dataset_key(), "within_composition") %||%
                    within_composition_choices(tab_object())
            })

            output[[paste0(prefix, "_distribution_group_by_ui")]] <- renderUI({
                choices <- tab_group_choices()

                selectInput(
                    inputId = paste0(prefix, "_distribution_group_by"),
                    label = "Color cells by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_distribution_group_by")]],
                        choices,
                        default = "Rank_1st"
                    )
                )
            })

            output[[paste0(prefix, "_distribution_split_by_ui")]] <- renderUI({
                choices <- tab_distribution_split_choices()

                selectInput(
                    inputId = paste0(prefix, "_distribution_split_by"),
                    label = "Split UMAP by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_distribution_split_by")]],
                        choices,
                        default = "none"
                    )
                )
            })

            output[[paste0(prefix, "_composition_by_ui")]] <- renderUI({
                choices <- tab_composition_choices()

                if (!length(choices)) {
                    return(NULL)
                }

                selectInput(
                    inputId = paste0(prefix, "_composition_by"),
                    label = "Show percentages by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_composition_by")]],
                        choices,
                        default = unname(choices[[1]])
                    )
                )
            })

            output[[paste0(prefix, "_split_by_ui")]] <- renderUI({
                choices <- tab_feature_split_choices()

                selectInput(
                    inputId = paste0(prefix, "_split_by"),
                    label = "Split feature UMAP by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_split_by")]],
                        choices,
                        default = "none"
                    )
                )
            })

            output[[paste0(prefix, "_notice_ui")]] <- renderUI({
                genes <- tab_local_applied_genes()

                if (!length(genes)) {
                    return(
                        div(
                            class = "alert-stack",
                            notice_card(
                                title = paste("No local", tab_label, "genes applied yet"),
                                body = sprintf(
                                    "Build a native %s gene panel to populate the expression plots in this tab. Cluster structure and composition remain available without a local panel.",
                                    tab_label
                                ),
                                tone = "info"
                            )
                        )
                    )
                }

                resolution <- tab_resolution()
                cards <- list()

                if (length(resolution$missing_features)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Genes missing from the selected atlas features",
                            body = paste(
                                compact_display_gene_list(species_key, resolution$missing_features, limit = 8),
                                "These genes are in the local panel, but matching features were not found in the current atlas object."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (isTRUE(resolution$feature_limit_exceeded)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Gene panel trimmed for rendering",
                            body = sprintf(
                                "This panel resolves to %s plotted features, so the app is showing the first %s to keep the browser responsive. Narrow the gene list or split it into smaller batches to inspect the omitted genes.",
                                format(resolution$feature_count_before_limit, big.mark = ","),
                                format(resolution$feature_limit, big.mark = ",")
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (!length(cards)) {
                    return(NULL)
                }

                div(class = "alert-stack", tagList(cards))
            })

            tab_distribution_group_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_distribution_group_by")]],
                    tab_group_choices(),
                    default = "Rank_1st"
                )
            })

            tab_distribution_split_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_distribution_split_by")]],
                    tab_distribution_split_choices(),
                    default = "none"
                )
            })

            tab_composition_by <- reactive({
                choices <- tab_composition_choices()
                default_choice <- if (length(choices)) unname(choices[[1]]) else NULL

                resolve_choice(
                    input[[paste0(prefix, "_composition_by")]],
                    choices,
                    default = default_choice
                )
            })

            tab_composition_cluster_by <- reactive({
                available_clusters <- clustering_columns[clustering_columns %in% unname(tab_group_choices())]
                preferred_cluster <- tab_distribution_group_by()

                if (length(preferred_cluster) && preferred_cluster %in% available_clusters) {
                    return(preferred_cluster)
                }

                if (length(available_clusters)) {
                    return(available_clusters[[1]])
                }

                NA_character_
            })

            tab_markers_full <- reactive({
                read_cluster_markers_cache(tab_dataset_key(), top_n = FALSE)
            })

            tab_marker_source <- reactive({
                marker_tbl <- tab_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NA_character_)
                }

                available_sources <- unique(as.character(marker_tbl$cluster_source))
                preferred_source <- tab_distribution_group_by()
                rank_label_columns <- paste0(setdiff(distribution_cluster_columns, "cluster_label"), "_label")
                if (length(preferred_source) && preferred_source %in% rank_label_columns) {
                    preferred_source <- sub("_label$", "", preferred_source)
                }
                if (!(length(preferred_source) && preferred_source %in% distribution_cluster_columns)) {
                    preferred_source <- NULL
                }

                select_marker_cluster_source(
                    available_sources = available_sources,
                    preferred_source = preferred_source,
                    default_source = "Rank_1st"
                )
            })

            tab_marker_lookup <- reactive({
                get_cached_cluster_lookup(tab_dataset_key(), tab_marker_source()) %||%
                    cluster_label_lookup(tab_object(), tab_marker_source())
            })

            tab_marker_choices <- reactive({
                marker_tbl <- tab_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(character(0))
                }

                marker_source <- tab_marker_source()

                marker_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source)

                available_clusters <- unique(as.character(marker_tbl$cluster))
                lookup_tbl <- tab_marker_lookup() %>%
                    filter(cluster %in% available_clusters)

                if (!nrow(lookup_tbl)) {
                    return(setNames(cluster_value_levels(available_clusters), cluster_value_levels(available_clusters)))
                }

                setNames(lookup_tbl$cluster, lookup_tbl$choice_label)
            })

            output[[paste0(prefix, "_marker_cluster_ui")]] <- renderUI({
                choices <- tab_marker_choices()

                if (!length(choices)) {
                    return(
                        tags$p(
                            class = "marker-status-hint",
                            tagList(
                                "Marker cache missing for this dataset. Run ",
                                tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                                "."
                            )
                        )
                    )
                }

                selectInput(
                    inputId = paste0(prefix, "_marker_cluster"),
                    label = "Cluster",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_marker_cluster")]],
                        choices,
                        default = unname(choices[[1]])
                    )
                )
            })

            observeEvent(tab_marker_choices(), {
                pending_cluster <- pending_url_cluster()

                if (is.null(pending_cluster) || !nzchar(pending_cluster) || !identical(input$main_tabs %||% "overview", species_key)) {
                    return()
                }

                choices <- tab_marker_choices()

                if (!length(choices)) {
                    return()
                }

                selected_cluster <- resolve_choice(
                    pending_cluster,
                    choices,
                    default = unname(choices[[1]])
                )

                updateSelectInput(
                    session,
                    inputId = paste0(prefix, "_marker_cluster"),
                    selected = selected_cluster
                )
                pending_url_cluster(NULL)
            }, ignoreInit = FALSE)

            tab_marker_cluster <- reactive({
                choices <- tab_marker_choices()
                default_choice <- if (length(choices)) unname(choices[[1]]) else NULL

                resolve_choice(
                    input[[paste0(prefix, "_marker_cluster")]],
                    choices,
                    default = default_choice
                )
            })

            tab_marker_top_n <- reactive({
                marker_n <- suppressWarnings(as.integer(input[[paste0(prefix, "_marker_top_n")]] %||% 10L))

                if (is.na(marker_n) || marker_n < 1L) {
                    marker_n <- 10L
                }

                min(marker_n, 25L)
            })

            tab_marker_table_raw <- reactive({
                marker_tbl <- tab_markers_full()

                validate(
                    need(
                        !is.null(marker_tbl),
                        tagList(
                            "Cluster markers are not cached for this dataset. Run ",
                            tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                            "."
                        )
                    ),
                    need(nrow(marker_tbl) > 0, "No cluster markers are available for this dataset.")
                )

                cluster_id <- tab_marker_cluster()
                marker_source <- tab_marker_source()
                filtered_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source) %>%
                    filter(cluster == !!cluster_id) %>%
                    arrange(p_val_adj, desc(avg_log2FC), desc(pct.1), gene)

                validate(need(nrow(filtered_tbl) > 0, "No marker rows are available for the selected cluster."))

                filtered_tbl
            })

            output[[paste0(prefix, "_markers_status_ui")]] <- renderUI({
                marker_tbl <- tab_markers_full()
                busy_message <- marker_job_messages()[[prefix]]

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NULL)
                }

                if (!is.null(busy_message) && nzchar(busy_message)) {
                    return(tags$p(class = "marker-status-hint is-working", busy_message))
                }

                cluster_choices <- tab_marker_choices()
                marker_source <- tab_marker_source()
                active_group_by <- tab_distribution_group_by()
                current_cluster_label <- names(cluster_choices)[match(tab_marker_cluster(), unname(cluster_choices))]
                current_cluster_label <- current_cluster_label[!is.na(current_cluster_label) & nzchar(current_cluster_label)][1] %||% tab_marker_cluster()
                source_sentence <- if (marker_cluster_source_matches(marker_source, active_group_by)) {
                    paste0("using ", marker_cluster_source_label(marker_source), ".")
                } else {
                    paste0(
                        "using ",
                        marker_cluster_source_label(marker_source),
                        " because ",
                        metadata_column_label(active_group_by),
                        " is not a cached clustering solution."
                    )
                }

                tags$p(
                    class = "marker-status-hint",
                    paste(
                        "Showing markers for",
                        current_cluster_label,
                        source_sentence,
                        "Top markers will be staged in this tab's Gene expression panel until you click Generate the expression plots."
                    )
                )
            })

            output[[paste0(prefix, "_markers_table")]] <- DT::renderDT({
                marker_tbl <- tab_marker_table_raw()

                display_tbl <- marker_tbl %>%
                    transmute(
                        gene = display_gene_labels(species_key, gene),
                        avg_log2FC = avg_log2FC,
                        pct.1 = pct.1,
                        pct.2 = pct.2,
                        p_val_adj = p_val_adj
                    )

                DT::datatable(
                    display_tbl,
                    rownames = FALSE,
                    escape = TRUE,
                    selection = "none",
                    class = "stripe hover order-column compact",
                    options = list(
                        dom = "tip",
                        pageLength = 10,
                        autoWidth = TRUE,
                        order = list(list(1, "desc"))
                    )
                ) %>%
                    DT::formatRound(columns = c("avg_log2FC", "pct.1", "pct.2"), digits = 3) %>%
                    DT::formatSignif(columns = "p_val_adj", digits = 3)
            })

            observeEvent(input[[paste0(prefix, "_add_markers")]], {
                button_id <- paste0(prefix, "_add_markers")
                marker_tbl <- tab_marker_table_raw() %>%
                    slice_head(n = tab_marker_top_n())

                set_marker_job_message(
                    prefix,
                    sprintf(
                        "Adding the top %d markers from %s into the local %s panel...",
                        nrow(marker_tbl),
                        tab_marker_cluster(),
                        tab_label
                    )
                )
                session$sendCustomMessage("atlas_button_busy", list(
                    button_id = button_id,
                    label = "Adding..."
                ))
                on.exit(
                    session$onFlushed(function() {
                        session$sendCustomMessage("atlas_button_busy", list(
                            button_id = button_id,
                            busy = FALSE
                        ))
                    }, once = TRUE),
                    add = TRUE
                )
                choice_bundle <- tab_local_choice_bundle()
                valid_local_ids <- choice_bundle$feature_ids %||% character(0)
                genes_to_add <- unique(marker_tbl$gene[marker_tbl$gene %in% valid_local_ids])

                if (!length(genes_to_add)) {
                    showNotification(
                        sprintf(
                            "No top markers from this cluster could be added to the local %s panel.",
                            tab_label
                        ),
                        type = "warning",
                        duration = 8
                    )
                    return()
                }

                updated_selection <- unique(c(tab_local_staged_genes(), genes_to_add))

                update_local_selected_genes_input(
                    species_key = species_key,
                    choice_bundle = choice_bundle,
                    selected = updated_selection
                )

                skipped_n <- nrow(marker_tbl) - length(genes_to_add)
                showNotification(
                    paste0(
                        "Added ",
                        length(genes_to_add),
                        " marker gene(s) to the local ",
                        tab_label,
                        " selection. Click Generate the expression plots to refresh this tab",
                        if (skipped_n > 0) paste0("; ", skipped_n, " were not present in this atlas.") else "."
                    ),
                    type = "message",
                    duration = 6
                )
            }, ignoreInit = TRUE)

            output[[paste0("dl_", prefix, "_markers")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_cluster_", tab_marker_cluster(), "_markers.csv")
                },
                content = function(file) {
                    marker_tbl <- tab_marker_table_raw()
                    lookup_tbl <- tab_marker_lookup()
                    cluster_labels <- lookup_tbl$cluster_label[match(marker_tbl$cluster, lookup_tbl$cluster)]

                    export_tbl <- marker_tbl %>%
                        mutate(
                            cluster_label = ifelse(is.na(cluster_labels) | !nzchar(cluster_labels), cluster, cluster_labels),
                            gene_label = display_gene_labels(species_key, gene)
                        ) %>%
                        select(cluster, cluster_label, gene, gene_label, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
                        add_export_provenance_columns(
                            tab_label = tab_label,
                            integration_label = current_species_integration(species_key),
                            extra = list(atlas_export_type = "cluster_markers")
                        )

                    write.csv(export_tbl, file = file, row.names = FALSE, na = "")
                }
            )

            tab_group_by <- reactive({
                resolve_choice(
                    tab_distribution_group_by(),
                    tab_group_choices(),
                    default = "Rank_1st"
                )
            })

            tab_split_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_split_by")]],
                    tab_feature_split_choices(),
                    default = "none"
                )
            })

            output[[paste0(prefix, "_distribution_umap_plot_ui")]] <- renderUI({
                split_by <- tab_distribution_split_by()

                if (identical(split_by, "none")) {
                    return(
                        div(
                            class = "distribution-split-layout",
                            div(
                                class = "distribution-split-panel",
                                div(class = "distribution-view-title", "2D UMAP"),
                                spinning_plot_output(
                                    paste0(prefix, "_distribution_umap_plot"),
                                    proxy_height = "540px",
                                    shell_class = "umap-plot-shell"
                                )
                            ),
                            div(
                                class = "distribution-split-panel",
                                div(class = "distribution-view-title", "3D UMAP"),
                                spinning_plotly_output(
                                    paste0(prefix, "_distribution_umap3d_plot"),
                                    proxy_height = "540px",
                                    shell_class = "plotly-plot-shell"
                                )
                            )
                        )
                    )
                }

                spinning_plot_output(
                    paste0(prefix, "_distribution_umap_plot"),
                    proxy_height = "520px",
                    shell_class = umap_plot_shell_class(split_by)
                )
            })

            output[[paste0(prefix, "_composition_plot_ui")]] <- renderUI({
                spinning_plot_output(
                    paste0(prefix, "_composition_plot"),
                    proxy_height = "480px"
                )
            })

            output[[paste0(prefix, "_umap_plot_ui")]] <- renderUI({
                spinning_plot_output(
                    paste0(prefix, "_umap_plot"),
                    proxy_height = "520px",
                    shell_class = "umap-plot-shell"
                )
            })

            distribution_umap_plot_obj <- reactive({
                obj <- tab_object()
                group_by <- tab_distribution_group_by()
                split_by <- tab_distribution_split_by()
                pt_size <- as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1)
                split_enabled <- !identical(split_by, "none")
                plot_obj <- apply_metadata_display_order(obj, c(group_by, split_by))
                color_map <- distribution_color_map(plot_obj@meta.data[[group_by]], group_by)
                colors_use <- unname(color_map)
                show_cluster_labels <- is_cluster_distribution_group(group_by) && identical(split_by, "none")
                split_panels <- if (split_enabled) split_panel_count(plot_obj, split_by) else 1L
                split_columns <- if (identical(split_by, "none")) {
                    NULL
                } else {
                    min(4L, split_panels)
                }

                distribution_plot <- scCustomize::DimPlot_scCustom(
                    seurat_object = plot_obj,
                    colors_use = colors_use,
                    group.by = group_by,
                    split.by = if (split_enabled) split_by else NULL,
                    pt.size = pt_size,
                    label = show_cluster_labels,
                    repel = TRUE,
                    raster = TRUE,
                    num_columns = split_columns
                )

                distribution_plot <- distribution_plot &
                    app_plot_theme() &
                    theme(
                        legend.title = element_blank(),
                        legend.position = if (split_enabled || is_cluster_distribution_group(group_by)) "none" else "top",
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        plot.margin = margin(8, 14, 10, 10)
                    )

                if (split_enabled) {
                    distribution_plot &
                        labs(color = NULL) &
                        theme(
                            plot.title = element_text(
                                face = "bold",
                                colour = app_palette["text"],
                                size = 13,
                                hjust = 0.5
                            )
                        )
                } else {
                    distribution_plot &
                        labs(title = NULL, color = NULL)
                }
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_distribution_group_by(),
                tab_distribution_split_by(),
                as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1),
                cache = "app"
            )

            distribution_umap3d_plot_data <- reactive({
                validate(
                    need(
                        identical(tab_distribution_split_by(), "none"),
                        "3D UMAP is available only when the distribution view is not split."
                    )
                )

                group_by <- tab_distribution_group_by()
                obj <- apply_metadata_display_order(tab_object(), group_by)
                embedding <- get_within_umap3d(
                    species_key = species_key,
                    integration_method = current_species_integration(species_key)
                )
                cell_ids <- intersect(rownames(obj@meta.data), rownames(embedding))

                validate(
                    need(length(cell_ids) > 0, "No 3D UMAP coordinates are available for this atlas.")
                )

                distribution_df <- tibble(
                    cell_id = cell_ids,
                    group_value = as.character(obj@meta.data[cell_ids, group_by, drop = TRUE]),
                    umap3d_1 = embedding[cell_ids, 1],
                    umap3d_2 = embedding[cell_ids, 2],
                    umap3d_3 = embedding[cell_ids, 3]
                ) %>%
                    filter(
                        !is.na(group_value) & nzchar(group_value),
                        !is.na(umap3d_1),
                        !is.na(umap3d_2),
                        !is.na(umap3d_3)
                    )

                validate(
                    need(nrow(distribution_df) > 0, "No cells are available for the 3D UMAP view.")
                )

                distribution_df <- distribution_df %>%
                    mutate(group_value = order_metadata_values(group_value, group_by))

                distribution_df <- stratified_point_sample(
                    distribution_df,
                    group_col = "group_value",
                    max_points = 30000L
                )

                full_color_map <- distribution_color_map(obj@meta.data[[group_by]], group_by)
                present_levels <- names(full_color_map)[names(full_color_map) %in% as.character(distribution_df$group_value)]

                list(
                    data = distribution_df,
                    color_map = full_color_map[present_levels],
                    group_by = group_by
                )
            })

            composition_plot_obj <- reactive({
                obj <- tab_object()
                composition_by <- tab_composition_by()
                cluster_by <- tab_composition_cluster_by()
                plot_obj <- apply_metadata_display_order(obj, composition_by)

                validate(
                    need(!is.null(composition_by) && nzchar(composition_by), "No composition metadata are available for this atlas."),
                    need(!is.na(cluster_by) && nzchar(cluster_by), "No clustering metadata are available for this atlas.")
                )

                md <- plot_obj@meta.data
                composition_df <- tibble(
                    cluster = as.character(md[[cluster_by]]),
                    composition = as.character(md[[composition_by]])
                ) %>%
                    filter(
                        !is.na(cluster) & nzchar(cluster),
                        !is.na(composition) & nzchar(composition)
                    ) %>%
                    count(cluster, composition, name = "cell_count")

                validate(
                    need(nrow(composition_df) > 0, "No cluster composition data are available for the selected metadata.")
                )

                composition_df <- composition_df %>%
                    mutate(
                        cluster = factor(cluster, levels = cluster_value_levels(cluster)),
                        composition = order_metadata_values(composition, composition_by)
                    )

                fill_values <- composition_colors_use(composition_df$composition, composition_by)
                legend_rows <- max(1L, min(3L, ceiling(dplyr::n_distinct(composition_df$composition) / 8L)))

                ggplot(
                    composition_df,
                    aes(x = cluster, y = cell_count, fill = composition)
                ) +
                    geom_col(
                        position = "fill",
                        width = 0.82,
                        colour = "#FFFFFF",
                        linewidth = 0.18
                    ) +
                    scale_y_continuous(
                        labels = scales::label_percent(accuracy = 1),
                        expand = expansion(mult = c(0, 0.02))
                    ) +
                    {if (!is.null(fill_values)) scale_fill_manual(values = fill_values)} +
                    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
                    labs(
                        x = metadata_column_label(cluster_by),
                        y = "Percent of cells",
                        fill = NULL
                    ) +
                    app_plot_theme() +
                    theme(
                        panel.grid.major.x = element_blank(),
                        axis.text.x = element_text(size = 11),
                        legend.position = "top",
                        legend.text = element_text(size = 9),
                        plot.margin = margin(8, 12, 8, 10)
                    )
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_composition_by(),
                tab_composition_cluster_by(),
                cache = "app"
            )

            rank_data <- reactive({
                resolution <- expression_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                dot_data <- Seurat::DotPlot(
                    object = obj,
                    features = resolution$plot_features,
                    group.by = tab_group_by()
                )$data

                label_map <- resolution$label_map

                as_tibble(dot_data) %>%
                    transmute(
                        gene_id = features.plot,
                        Gene = unname(label_map[features.plot]),
                        Group = id,
                        `Pct. expressing` = round(pct.exp, 1),
                        `Scaled expression` = round(avg.exp.scaled, 2),
                        `Average expression` = round(avg.exp, 3)
                    ) %>%
                    group_by(Gene) %>%
                    slice_max(order_by = `Scaled expression`, n = 5, with_ties = FALSE) %>%
                    ungroup()
            })

	            umap_plot_obj <- reactive({
	                resolution <- expression_resolution()
	                obj <- tab_object()
	                source_species <- species_key
	                reference_group_by <- tab_composition_cluster_by()
	                split_by <- tab_split_by()
	                pt_size <- max(2.0, as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1) * 1.8)
	                split_enabled <- !identical(split_by, "none")
	                feature_panel_n <- length(resolution$plot_features) + 1L
	                plot_obj <- apply_metadata_display_order(obj, c(split_by, reference_group_by))
	                split_panels <- if (split_enabled) split_panel_count(plot_obj, split_by) else 1L
	                split_columns <- if (split_enabled) min(4L, split_panels) else NULL
	                feature_grid_cols <- within_feature_grid_cols(
	                    feature_n = feature_panel_n,
	                    split_by = split_by
	                )

	                validate(
	                    need(
	                        length(resolution$plot_features) > 0,
	                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
	                    ),
	                    need(
	                        !is.na(reference_group_by) && nzchar(reference_group_by),
	                        paste("No clustering metadata are available for", tab_label, ".")
	                    )
	                )

	                reference_color_map <- distribution_color_map(
	                    plot_obj@meta.data[[reference_group_by]],
	                    reference_group_by
	                )

	                reference_plot <- scCustomize::DimPlot_scCustom(
	                    seurat_object = plot_obj,
	                    colors_use = unname(reference_color_map),
	                    group.by = reference_group_by,
	                    split.by = if (split_enabled) split_by else NULL,
	                    pt.size = pt_size,
	                    label = identical(split_by, "none"),
	                    repel = TRUE,
	                    raster = TRUE,
	                    num_columns = split_columns
	                )

	                reference_plot <- reference_plot &
	                    labs(color = NULL) &
	                    app_plot_theme() &
	                    theme(
	                        legend.title = element_blank(),
	                        legend.position = "none",
	                        panel.grid = element_blank(),
	                        axis.title = element_blank(),
	                        axis.text = element_blank(),
	                        axis.ticks = element_blank(),
	                        axis.line = element_blank(),
	                        plot.margin = margin(4, 8, 6, 6)
	                    )

	                reference_plot <- if (split_enabled) {
	                    reference_plot &
	                        theme(
	                            plot.title = element_text(
	                                face = "bold",
	                                colour = app_palette["text"],
	                                size = 13,
	                                hjust = 0.5
	                            )
	                        )
	                } else {
	                    reference_plot &
	                        labs(title = NULL) &
	                        theme(plot.title = element_blank())
	                }

	                plot_list <- c(list(
	                    wrap_titled_plot(
	                        plot_obj = reference_plot,
	                        title = metadata_column_label(reference_group_by)
	                    )
	                ), lapply(resolution$plot_features, function(feature_id) {
	                    feature_plot_args <- list(
	                        seurat_object = plot_obj,
	                        features = feature_id,
                        split.by = if (split_enabled) split_by else NULL,
                        pt.size = pt_size,
                        order = TRUE,
                        raster = TRUE,
                        num_columns = split_columns,
                        label = FALSE,
                        combine = TRUE
                    )
                    if (isTRUE(input$colorblind_safe)) {
                        feature_plot_args$colors_use <- viridisLite::cividis(100, end = 0.95)
                    }
                    feature_plot <- do.call(scCustomize::FeaturePlot_scCustom, feature_plot_args)

                    feature_plot <- feature_plot &
                        labs(color = NULL) &
                        app_plot_theme() &
                        compact_feature_legend_guides() &
                        compact_feature_legend_theme() &
                        scale_x_continuous(expand = expansion(mult = 0.01)) &
                        scale_y_continuous(expand = expansion(mult = 0.01)) &
                        theme(
                            panel.grid = element_blank(),
                            axis.title = element_blank(),
                            axis.text = element_blank(),
                            axis.ticks = element_blank(),
                            axis.line = element_blank(),
                            plot.margin = margin(4, 8, 6, 6)
                        )

                    feature_plot <- if (split_enabled) {
                        feature_plot &
                            theme(
                                plot.title = element_text(
                                    face = "bold",
                                    colour = app_palette["text"],
                                    size = 13,
                                    hjust = 0.5
                                )
                            )
                    } else {
                        feature_plot &
                            labs(title = NULL) &
                            theme(plot.title = element_blank())
                    }

	                    wrap_titled_plot(
	                        plot_obj = feature_plot,
	                        title = format_within_feature_panel_title(
	                            title = unname(resolution$label_map[feature_id]),
	                            source_species = source_species,
	                            target_species = species_key
	                        )
	                    )
	                }))

	                wrap_plots(plotlist = plot_list, ncol = feature_grid_cols)
	            # Cache local plots against this tab's applied local genes so
	            # the within-species Generate button always invalidates them.
	            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_local_applied_genes(),
                tab_split_by(),
                tab_composition_cluster_by(),
                as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            violin_plot_obj <- reactive({
                resolution <- expression_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                label_map <- resolution$label_map
                n_groups <- dplyr::n_distinct(obj@meta.data[[tab_group_by()]])

                plot_list <- lapply(resolution$plot_features, function(feature_id) {
                    violin_plot <- scCustomize::VlnPlot_scCustom(
                        seurat_object = obj,
                        features = feature_id,
                        group.by = tab_group_by(),
                        colors_use = rep(unname(app_palette["warm"]), n_groups),
                        pt.size = 0.12,
                        num_columns = 1,
                        raster = FALSE
                    )

                    if (length(violin_plot$layers) >= 2) {
                        violin_plot$layers[[2]]$aes_params$colour <- unname(app_palette["green_dark"])
                        violin_plot$layers[[2]]$aes_params$alpha <- 0.42
                        violin_plot$layers[[2]]$aes_params$size <- 0.34
                    }

                    violin_plot <- violin_plot &
                        labs(
                            title = unname(label_map[feature_id]),
                            x = NULL,
                            y = "Normalized expression"
                        ) &
                        app_plot_theme() &
                        theme(
                            plot.title = element_text(
                                face = "bold",
                                colour = app_palette["text"],
                                size = 16,
                                hjust = 0
                            ),
                            legend.position = "none",
                            axis.text.x = element_text(angle = 35, hjust = 1),
                            panel.grid.major.x = element_blank(),
                            plot.margin = margin(8, 14, 12, 10)
                        )

                    violin_plot
                })

                if (length(plot_list) == 1) {
                    plot_list[[1]]
                } else {
                    wrap_plots(plotlist = plot_list, ncol = 1)
                }
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_local_applied_genes(),
                tab_group_by(),
                cache = "app"
            )

            heatmap_plot_obj <- reactive({
                resolution <- expression_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                build_expression_heatmap_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = tab_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_local_applied_genes(),
                tab_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            ridge_plot_obj <- reactive({
                resolution <- expression_resolution()
                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

                build_expression_ridge_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = tab_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_local_applied_genes(),
                tab_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

	            dot_plot_obj <- reactive({
	                resolution <- expression_resolution()
	                obj <- tab_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste("No mapped genes are available for", tab_label, "in the selected atlas.")
                    )
                )

	                label_map <- unname(resolution$label_map[resolution$plot_features])
	                cluster_feature_enabled <- length(resolution$plot_features) > 1

	                dot_plot_precheck <- tryCatch(
	                    scCustomize::DotPlot_scCustom(
	                        seurat_object = obj,
	                        features = resolution$plot_features,
	                        group.by = tab_group_by()
	                    )$data,
	                    error = function(e) NULL
	                )

	                if (!is.null(dot_plot_precheck) && "avg.exp.scaled" %in% colnames(dot_plot_precheck)) {
	                    has_nonfinite_feature <- any(
	                        !is.finite(dot_plot_precheck$avg.exp.scaled) &
	                            !is.na(dot_plot_precheck$features.plot)
	                    )

	                    if (isTRUE(has_nonfinite_feature)) {
	                        cluster_feature_enabled <- FALSE
	                    }
	                }

	                build_clustered_dot_plot <- function(cluster_feature) {
	                    scCustomize::Clustered_DotPlot(
	                        seurat_object = obj,
	                        features = resolution$plot_features,
	                        group.by = tab_group_by(),
	                        cluster_feature = cluster_feature,
	                        cluster_ident = TRUE,
	                        flip = FALSE,
	                        row_names_side = "left",
	                        column_names_side = "bottom",
	                        show_ident_colors = FALSE,
	                        show_ident_legend = FALSE,
	                        grid_color = NA,
	                        row_label_size = 11,
	                        column_label_size = 11,
	                        legend_position = "right",
	                        legend_label_size = 10,
	                        legend_title_size = 10,
	                        x_lab_rotate = TRUE,
	                        raster = FALSE,
	                        plot_km_elbow = FALSE
	                    )
	                }

	                plot_obj <- tryCatch(
	                    build_clustered_dot_plot(cluster_feature_enabled),
	                    error = function(e) {
	                        if (isTRUE(cluster_feature_enabled)) {
	                            build_clustered_dot_plot(FALSE)
	                        } else {
	                            stop(e)
	                        }
	                    }
	                )

	                clustered_dot_radius_mm <- 3.2
	                legend_dot_radius_mm <- 2
	                heatmap_obj <- plot_obj@ht_list[[1]]
	                percent_scale_max <- NA_real_

	                if (!is.null(heatmap_obj@matrix_param$cell_fun)) {
	                    cell_fun_env <- environment(heatmap_obj@matrix_param$cell_fun)
	                    percent_mat <- get("percent_mat", envir = cell_fun_env)
	                    percent_scale_max <- max(percent_mat, na.rm = TRUE)
	                    heatmap_obj@matrix_param$cell_fun <- eval(
	                        bquote(
	                            function(j, i, x, y, w, h, fill) {
	                                grid::grid.rect(
	                                    x = x,
	                                    y = y,
	                                    width = w,
	                                    height = h,
	                                    gp = grid::gpar(col = grid_color, fill = NA)
	                                )
	                                grid::grid.circle(
	                                    x = x,
	                                    y = y,
	                                    r = sqrt(percent_mat[i, j] / .(max(percent_scale_max, 1e-08))) * grid::unit(.(clustered_dot_radius_mm), "mm"),
	                                    gp = grid::gpar(fill = col_fun(exp_mat[i, j]), col = NA)
	                                )
	                            }
	                        ),
	                        envir = cell_fun_env
	                    )
	                }

	                if (!is.null(heatmap_obj@matrix_param$layer_fun)) {
	                    layer_fun_env <- environment(heatmap_obj@matrix_param$layer_fun)
	                    if (!is.finite(percent_scale_max) || percent_scale_max <= 0) {
	                        percent_scale_max <- max(get("percent_mat", envir = layer_fun_env), na.rm = TRUE)
	                    }
	                    heatmap_obj@matrix_param$layer_fun <- eval(
	                        bquote(
	                            function(j, i, x, y, w, h, fill) {
	                                grid::grid.rect(
	                                    x = x,
	                                    y = y,
	                                    width = w,
	                                    height = h,
	                                    gp = grid::gpar(col = grid_color, fill = NA)
	                                )
	                                grid::grid.circle(
	                                    x = x,
	                                    y = y,
	                                    r = sqrt(ComplexHeatmap::pindex(percent_mat, i, j) / .(max(percent_scale_max, 1e-08))) * grid::unit(.(clustered_dot_radius_mm), "mm"),
	                                    gp = grid::gpar(
	                                        fill = col_fun(ComplexHeatmap::pindex(exp_mat, i, j)),
	                                        col = NA
	                                    )
	                                )
	                            }
	                        ),
	                        envir = layer_fun_env
	                    )
	                }

	                if (!is.finite(percent_scale_max) || percent_scale_max <= 0) {
	                    percent_scale_max <- 1
	                }

	                legend_breaks <- pretty(c(0, percent_scale_max), n = 4)
	                legend_breaks <- legend_breaks[legend_breaks > 0 & legend_breaks < percent_scale_max]
	                legend_breaks <- unique(c(legend_breaks, percent_scale_max))
	                legend_breaks <- legend_breaks[legend_breaks > 0]

	                if (!length(legend_breaks)) {
	                    legend_breaks <- percent_scale_max
	                }

	                legend_digits <- if (percent_scale_max < 5) {
	                    2L
	                } else if (percent_scale_max < 20) {
	                    1L
	                } else {
	                    0L
	                }

	                legend_labels <- sub(
	                    "\\.?0+$",
	                    "",
	                    formatC(legend_breaks, format = "f", digits = legend_digits)
	                )

	                size_legend <- ComplexHeatmap::Legend(
	                    title = "Percent Expressing",
	                    at = legend_labels,
	                    graphics = lapply(legend_breaks, function(value) {
	                        force(value)
	                        function(x, y, w, h) {
	                            grid::grid.circle(
	                                x = x,
	                                y = y,
	                                r = sqrt(value / percent_scale_max) * grid::unit(legend_dot_radius_mm, "mm"),
	                                gp = grid::gpar(fill = "black", col = NA)
	                            )
	                        }
	                    }),
	                    labels_gp = grid::gpar(fontsize = 10),
	                    title_gp = grid::gpar(fontsize = 10, fontface = "bold")
	                )

	                row_feature_ids <- rownames(heatmap_obj@matrix)
	                if (is.null(row_feature_ids) || !length(row_feature_ids)) {
	                    row_feature_ids <- resolution$plot_features
	                }

	                dot_plot_row_labels <- unname(resolution$label_map[row_feature_ids])
	                missing_label_idx <- which(is.na(dot_plot_row_labels) | !nzchar(dot_plot_row_labels))
	                if (length(missing_label_idx)) {
	                    dot_plot_row_labels[missing_label_idx] <- row_feature_ids[missing_label_idx]
	                }

	                plot_obj@ht_list[[1]] <- set_heatmap_row_labels(
	                    heatmap_obj,
	                    dot_plot_row_labels
	                )

	                if (length(plot_obj@heatmap_legend_param$list)) {
	                    plot_obj@heatmap_legend_param$list[[1]] <- size_legend
	                } else {
	                    plot_obj@heatmap_legend_param$list <- list(size_legend)
	                }

	                plot_obj
	            }) %>% bindCache(
                species_key,
                current_species_integration(species_key),
                tab_local_applied_genes(),
                tab_group_by(),
                cache = "app"
            )

            draw_clustered_dot_plot <- function(plot_obj) {
                ComplexHeatmap::draw(
                    plot_obj,
                    heatmap_legend_side = "right",
                    annotation_legend_side = "right",
                    merge_legends = TRUE
                )
            }

            clustered_dot_plot_height <- function() {
                feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                group_n <- tryCatch(dplyr::n_distinct(tab_object()@meta.data[[tab_group_by()]]), error = function(e) 0L)

                max(
                    460,
                    170 + feature_n * 42 + group_n * 11
                )
            }

            output[[paste0(prefix, "_umap_plot")]] <- renderPlot(
                {
                    validate_local_expression_ready()
                    umap_plot_obj()
                },
	                height = function() {
	                    feature_n <- tryCatch(length(expression_resolution()$plot_features) + 1L, error = function(e) 0L)
	                    split_by <- tryCatch(tab_split_by(), error = function(e) "none")
	                    feature_cols <- if (identical(split_by, "none")) {
	                        within_feature_grid_cols(feature_n = feature_n, split_by = split_by)
                    } else {
                        4L
                    }

                    if (identical(split_by, "none")) {
                        feature_umap_height_px(
                            feature_n = feature_n,
                            feature_cols = feature_cols,
                            split_by = split_by
                        )
                    } else {
                        panels_per_gene <- tryCatch(
                            split_panel_count(tab_object(), split_by),
                            error = function(e) 1L
                        )
                        feature_umap_height_px(
                            feature_n = feature_n,
                            feature_cols = feature_cols,
                            split_by = split_by,
                            panels_per_gene = panels_per_gene
                        )
                    }
                },
                res = 110
            )

            output[[paste0(prefix, "_distribution_umap_plot")]] <- renderPlot(
                {
                    distribution_umap_plot_obj()
                },
                height = function() {
                    split_by <- tryCatch(tab_distribution_split_by(), error = function(e) "none")

                    if (identical(split_by, "none")) {
                        return(540)
                    }

                    panel_n <- tryCatch(
                        split_panel_count(tab_object(), split_by),
                        error = function(e) 1L
                    )
                    split_cols <- min(4L, max(1L, panel_n))
                    split_rows <- ceiling(panel_n / split_cols)

                    max(520, split_rows * 320)
                },
                res = 110
            )

            output[[paste0(prefix, "_distribution_umap3d_plot")]] <- plotly::renderPlotly({
                plot_data <- distribution_umap3d_plot_data()
                df <- plot_data$data
                color_map <- plot_data$color_map
                group_by <- plot_data$group_by
                marker_size <- max(1.8, as.numeric(input[[paste0(prefix, "_distribution_pt_size")]] %||% 1.1) * 2.2)
                group_levels <- names(color_map)

                p <- plotly::plot_ly(
                    type = "scatter3d",
                    mode = "markers"
                )

                for (group_level in group_levels) {
                    group_df <- df %>%
                        filter(as.character(group_value) == !!group_level)

                    if (!nrow(group_df)) {
                        next
                    }

                    p <- p %>%
                        plotly::add_trace(
                            data = group_df,
                            x = ~umap3d_1,
                            y = ~umap3d_2,
                            z = ~umap3d_3,
                            name = group_level,
                            marker = list(
                                size = marker_size,
                                color = unname(color_map[group_level]),
                                opacity = 0.8
                            ),
                            text = ~paste0(
                                metadata_column_label(group_by), ": ", group_value,
                                "<br>Cell: ", cell_id
                            ),
                            hovertemplate = "%{text}<extra></extra>"
                        )
                }

                p %>%
                    plotly::layout(
                        margin = list(l = 0, r = 0, b = 0, t = 10),
                        legend = list(
                            orientation = "h",
                            x = 0,
                            y = 1.08,
                            font = list(
                                size = 14,
                                color = unname(app_palette["text"])
                            ),
                            itemsizing = "constant",
                            itemwidth = 40,
                            bgcolor = "rgba(255,255,255,0.82)",
                            bordercolor = "rgba(201,214,196,0.92)",
                            borderwidth = 1
                        ),
                        scene = list(
                            aspectmode = "data",
                            xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            camera = list(
                                eye = list(x = 1.22, y = 1.06, z = 0.75)
                            )
                        )
                    ) %>%
                    plotly::config(displaylogo = FALSE)
            })

            output[[paste0(prefix, "_composition_plot")]] <- renderPlot(
                {
                    composition_plot_obj()
                },
                height = function() {
                    composition_by <- tryCatch(tab_composition_by(), error = function(e) NULL)
                    obj <- tryCatch(tab_object(), error = function(e) NULL)

                    if (is.null(obj) || is.null(composition_by) || !nzchar(composition_by)) {
                        return(460)
                    }

                    level_count <- tryCatch({
                        values <- obj@meta.data[[composition_by]]
                        values <- as.character(values)
                        length(unique(values[!is.na(values) & nzchar(values)]))
                    }, error = function(e) 1L)

                    legend_rows <- max(1L, min(3L, ceiling(level_count / 8L)))
                    max(460, 380 + legend_rows * 44)
                },
                res = 110
            )

            output[[paste0(prefix, "_violin_plot")]] <- renderPlot(
                {
                    validate_local_expression_ready()
                    violin_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                    max(320, 250 * feature_n)
                },
                res = 110
            )

            output[[paste0(prefix, "_heatmap_plot")]] <- renderPlot(
                {
                    validate_local_expression_ready()
                    heatmap_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                    max(420, 110 + feature_n * 34)
                },
                res = 110
            )

            output[[paste0(prefix, "_ridge_plot")]] <- renderPlot(
                {
                    validate_local_expression_ready()
                    ridge_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(expression_resolution()$plot_features), error = function(e) 0L)
                    group_n <- tryCatch(dplyr::n_distinct(tab_object()@meta.data[[tab_group_by()]]), error = function(e) 0L)
                    expression_ridge_height_px(feature_n = feature_n, group_n = group_n)
                },
                res = 110
            )

            output[[paste0(prefix, "_dot_plot")]] <- renderPlot(
                {
                    validate_local_expression_ready()
                    draw_clustered_dot_plot(dot_plot_obj())
                },
                height = function() {
                    clustered_dot_plot_height()
                },
                res = 110
            )

	            output[[paste0("dl_", prefix, "_umap")]] <- downloadHandler(
		                filename = function() {
		                    paste0(prefix, "_umap.", get_ext())
		                },
	                content = function(file) {
	                    resolution <- expression_resolution()
	                    feature_n <- length(resolution$plot_features) + 1L
	                    split_by <- tab_split_by()
	                    feature_cols <- if (identical(split_by, "none")) {
	                        within_feature_grid_cols(feature_n = feature_n, split_by = split_by)
                    } else {
                        4L
                    }

                    plot_height <- feature_umap_height_inches(
                        feature_n = feature_n,
                        feature_cols = feature_cols,
                        split_by = split_by,
                        panels_per_gene = if (identical(split_by, "none")) {
                            1L
                        } else {
                            split_panel_count(tab_object(), split_by)
                        }
                    )

                    save_ggplot(
                        file = file,
                        plot_obj = umap_plot_obj(),
                        width = 15,
                        height = plot_height,
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "umap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_distribution_umap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_distribution_umap.", get_ext())
                },
                content = function(file) {
                    split_by <- tab_distribution_split_by()

                    plot_height <- if (identical(split_by, "none")) {
                        7
                    } else {
                        panel_n <- split_panel_count(tab_object(), split_by)
                        split_cols <- min(4L, max(1L, panel_n))
                        split_rows <- ceiling(panel_n / split_cols)
                        max(7, split_rows * 3.4)
                    }

                    save_ggplot(
                        file = file,
                        plot_obj = distribution_umap_plot_obj(),
                        width = 14,
                        height = plot_height,
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "distribution_umap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_violin")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_violins.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(expression_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = violin_plot_obj(),
                        width = 10,
                        height = max(6, feature_n * 3.2),
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "violin")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_heatmap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_heatmap.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(expression_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = heatmap_plot_obj(),
                        width = 13,
                        height = max(5.5, 1.8 + feature_n * 0.55),
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "heatmap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_ridge")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_ridgeplot.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(expression_resolution()$plot_features)
                    group_n <- dplyr::n_distinct(tab_object()@meta.data[[tab_group_by()]])

                    save_ggplot(
                        file = file,
                        plot_obj = ridge_plot_obj(),
                        width = 10,
                        height = max(6, expression_ridge_height_px(feature_n, group_n) / 95),
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "ridgeplot")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_dot")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_dotplot.", get_ext())
                },
                content = function(file) {
                    plot_obj <- dot_plot_obj()
                    ext <- get_ext()
                    plot_height <- clustered_dot_plot_height() / 95
                    preset_cfg <- figure_preset_config(figure_preset())
                    export_caption <- export_plot_caption(
                        tab_label = tab_label,
                        integration_label = current_species_integration(species_key),
                        extra = list(atlas_export_type = "dotplot")
                    )

                    if (identical(ext, "png")) {
                        png(
                            filename = file,
                            width = 1900 * preset_cfg$width_scale,
                            height = max(1200, clustered_dot_plot_height() * 2) * preset_cfg$height_scale,
                            res = 220
                        )
                    } else if (identical(ext, "svg")) {
                        svglite::svglite(
                            file = file,
                            width = 14 * preset_cfg$width_scale,
                            height = plot_height * preset_cfg$height_scale
                        )
                    } else {
                        pdf(
                            file = file,
                            width = 14 * preset_cfg$width_scale,
                            height = plot_height * preset_cfg$height_scale
                        )
                    }

                    draw_clustered_dot_plot(plot_obj)
                    grid::grid.text(
                        export_caption,
                        x = grid::unit(0.01, "npc"),
                        y = grid::unit(0.01, "npc"),
                        just = c("left", "bottom"),
                        gp = grid::gpar(fontsize = 8.5, col = unname(app_palette["muted"]))
                    )
                    grDevices::dev.off()
                }
            )
        })
    }

    walk(within_species_keys, register_species_tab)

    register_cross_tab <- function(cross_key) {
        local({
            integration_cfg <- cross_integration_registry[[cross_key]]
            prefix <- paste0("cross_", cross_key)
            cross_tab_is_active <- reactive({
                identical(input$main_tabs %||% "overview", prefix)
            })

            cross_object <- reactive({
                req(cross_tab_is_active())
                get_cross_object(cross_key)
            })

            cross_resolution <- reactive({
                req(cross_tab_is_active())
                enforce_plot_feature_limit(resolve_cross_integration_mapping(
                    source_species = current_comparison_source_species(),
                    source_genes = selected_source_genes(),
                    cross_key = cross_key
                ))
            }) %>% bindCache(
                current_comparison_source_species(),
                selected_source_genes(),
                cross_key,
                atlas_plot_feature_limit,
                cache = "app"
            )

            output[[paste0(prefix, "_comparison_panel_ui")]] <- renderUI({
                if (!isTRUE(cross_tab_is_active())) {
                    return(NULL)
                }

                tagList(
                    div(
                        class = "subsection-header",
                        h3("Cross-species comparison panel"),
                        p("Choose the source species for this comparison, stage genes from that species, and apply the shared panel used throughout this integration tab.")
                    ),
                    fluidRow(
                        column(
                            width = 3,
                            div(
                                class = "option-group source-species-picker",
                                shinyWidgets::pickerInput(
                                    inputId = "source_species",
                                    label = "Comparison source species",
                                    choices = species_choices,
                                    selected = current_comparison_source_species(),
                                    multiple = FALSE
                                )
                            )
                        ),
                        column(
                            width = 6,
                            div(
                                class = "option-group",
                                selectizeInput(
                                    inputId = "selected_genes",
                                    label = "Comparison panel genes",
                                    choices = NULL,
                                    multiple = TRUE,
                                    options = gene_panel_selectize_options(enable_bulk = TRUE)
                                ),
                                div(
                                    class = "gene-action-row",
                                    actionButton(
                                        inputId = "apply_gene_selection",
                                        label = "Apply comparison panel",
                                        icon = icon("play"),
                                        class = "btn btn-default apply-selection-btn"
                                    ),
                                    actionButton(
                                        inputId = "open_gene_import",
                                        label = "Import list...",
                                        icon = icon("file-import"),
                                        class = "btn btn-default gene-import-btn"
                                    )
                                ),
                                div(
                                    class = "selection-meta",
                                    `aria-live` = "polite",
                                    `aria-atomic` = "true",
                                    textOutput("gene_selection_status")
                                )
                            )
                        ),
                        column(
                            width = 3,
                            div(
                                class = "option-group",
                                tags$p(
                                    class = "marker-status-hint",
                                    "The selected genes remain shared between Camex and SATURN so you can switch integrations without rebuilding the cross-species comparison."
                                )
                            )
                        )
                    )
                )
            })

            observe({
                if (!isTRUE(cross_tab_is_active())) {
                    return()
                }

                update_selected_genes_input(
                    choice_bundle = build_gene_choices(
                        current_comparison_source_species(),
                        source_integration()
                    ),
                    selected = staged_source_genes()
                )
            })

            cross_group_choices_cached <- reactive({
                get_cached_ui_choices(cross_key, "cross_group") %||%
                    cross_group_choices(cross_object())
            })

            cross_distribution_group_choices_cached <- reactive({
                get_cached_ui_choices(cross_key, "cross_distribution_group") %||%
                    cross_distribution_group_choices(cross_object(), cross_key)
            })

            cross_composition_choices_cached <- reactive({
                get_cached_ui_choices(cross_key, "cross_composition") %||%
                    cross_composition_choices(cross_object())
            })

            output[[paste0(prefix, "_group_by_ui")]] <- renderUI({
                choices <- cross_group_choices_cached()

                selectInput(
                    inputId = paste0(prefix, "_group_by"),
                    label = "Summarize dot plot by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_group_by")]],
                        choices,
                        default = integration_cfg$default_group_by
                    )
                )
            })

            cross_group_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_group_by")]],
                    cross_group_choices_cached(),
                    default = integration_cfg$default_group_by
                )
            })

            cross_dist_group_by <- reactive({
                resolve_choice(
                    input[[paste0(prefix, "_dist_group_by")]],
                    cross_distribution_group_choices_cached(),
                    default = integration_cfg$default_group_by
                )
            })

            output[[paste0(prefix, "_dist_group_by_ui")]] <- renderUI({
                choices <- cross_distribution_group_choices_cached()
                selectInput(
                    inputId = paste0(prefix, "_dist_group_by"),
                    label = "Color cells by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_dist_group_by")]],
                        choices,
                        default = integration_cfg$default_group_by
                    )
                )
            })

            output[[paste0(prefix, "_composition_by_ui")]] <- renderUI({
                choices <- cross_composition_choices_cached()
                if (!length(choices)) return(NULL)
                selectInput(
                    inputId = paste0(prefix, "_composition_by"),
                    label = "Show percentages by",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_composition_by")]],
                        choices,
                        default = unname(choices[[1]])
                    )
                )
            })

            output[[paste0(prefix, "_dist_umap_plot_ui")]] <- renderUI({
                has_3d <- !is.null(get_cross_umap3d(cross_key))

                if (has_3d) {
                    div(
                        class = "distribution-split-layout",
                        div(
                            class = "distribution-split-panel",
                            div(class = "distribution-view-title", "2D UMAP"),
                            spinning_plot_output(
                                paste0(prefix, "_dist_umap_plot"),
                                proxy_height = "540px",
                                shell_class = "umap-plot-shell"
                            )
                        ),
                        div(
                            class = "distribution-split-panel",
                            div(class = "distribution-view-title", "3D UMAP"),
                            spinning_plotly_output(
                                paste0(prefix, "_dist_umap3d_plot"),
                                proxy_height = "540px",
                                shell_class = "plotly-plot-shell"
                            )
                        )
                    )
                } else {
                    spinning_plot_output(
                        paste0(prefix, "_dist_umap_plot"),
                        proxy_height = "540px",
                        shell_class = "umap-plot-shell"
                    )
                }
            })

            output[[paste0(prefix, "_composition_plot_ui")]] <- renderUI({
                spinning_plot_output(
                    paste0(prefix, "_composition_plot"),
                    proxy_height = "480px"
                )
            })

            output[[paste0(prefix, "_umap_plot_ui")]] <- renderUI({
                spinning_plot_output(
                    paste0(prefix, "_umap_plot"),
                    proxy_height = "520px",
                    shell_class = "umap-plot-shell"
                )
            })

            output[[paste0(prefix, "_notice_ui")]] <- renderUI({
                genes <- selected_source_genes()
                resolution <- cross_resolution()
                feature_mode <- integration_cfg$feature_mode
                source_species <- current_comparison_source_species()
                render_summary <- cross_render_load_summary(resolution, source_species, cross_key)

                intro_card <- if (identical(feature_mode, "medicago_space")) {
                    notice_card(
                        title = "Shared Medicago-space features",
                        body = paste(
                            cross_integration_label(cross_key),
                            "stores a shared feature space represented with Medicago identifiers. Soybean and Lotus selections are projected to Medicago orthologs in this tab, and neighborhood overlap in the integrated embedding should be treated as comparative structure rather than proof of conserved cell states."
                        ),
                        tone = "info"
                    )
                } else {
                    notice_card(
                        title = paste(cross_integration_label(cross_key), "ortholog feature space"),
                        body = paste(
                            cross_integration_label(cross_key),
                            "stores species-prefixed features. Source genes are resolved to ortholog features from Medicago, Glycine, and Lotus before the integrated plots are drawn, and the shared embedding should be interpreted as a comparative view rather than a one-to-one cell-state map."
                        ),
                        tone = "info"
                    )
                }

                cards <- list(intro_card)

                if (!length(genes)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "No source genes selected",
                            body = paste(
                                "Add one or more source-species genes to populate the",
                                cross_integration_label(cross_key),
                                "plots."
                            ),
                            tone = "info"
                        )
                    ))
                }

                if (length(resolution$no_orthogroup)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Selected genes without orthogroups",
                            body = paste(
                                compact_display_gene_list(current_comparison_source_species(), resolution$no_orthogroup, limit = 8),
                                "These genes are not represented in the current orthogroup table, so the cross-species panels cannot speak to them."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (length(resolution$no_target_members)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = if (identical(feature_mode, "medicago_space")) {
                                "Orthogroups without Medicago members"
                            } else {
                                "Orthogroups without mapped members"
                            },
                            body = if (identical(feature_mode, "medicago_space")) {
                                compact_display_gene_list(current_comparison_source_species(), resolution$no_target_members, limit = 8)
                            } else {
                                paste(
                                    "No ortholog members from Medicago, Glycine, or Lotus were found in the mapped orthogroups for:",
                                    compact_display_gene_list(current_comparison_source_species(), resolution$no_target_members, limit = 8),
                                    "Treat this as an orthogroup-content limit, not as evidence that the biology is absent in the other species."
                                )
                            },
                            tone = "warning"
                        )
                    ))
                }

                if (length(resolution$missing_features)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = paste("Mapped orthologs missing from the", cross_integration_label(cross_key), "feature set"),
                            body = paste(
                                compact_display_gene_list(current_comparison_source_species(), resolution$missing_features, limit = 8),
                                "Orthologs were resolved but the integration feature set does not contain them. Treat this as feature-space coverage loss rather than proof of no expression."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (nrow(resolution$multiplicity)) {
                    multiplicity_text <- resolution$multiplicity %>%
                        mutate(
                            label = paste0(
                                display_gene_labels(source_species, source_gene, include_gene_id_with_common = FALSE),
                                " (", mapped_gene_count, " mapped features)"
                            )
                        ) %>%
                        pull(label)

                    cards <- append(cards, list(
                        notice_card(
                            title = paste("One-to-many mappings in", cross_integration_label(cross_key)),
                            body = paste(
                                compact_gene_list(multiplicity_text, limit = 6),
                                "These are ambiguous orthogroup mappings. Each mapped gene is plotted separately, and none should be treated as the definitive conserved counterpart without outside support."
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (isTRUE(resolution$feature_limit_exceeded)) {
                    cards <- append(cards, list(
                        notice_card(
                            title = "Mapped feature panel trimmed for rendering",
                            body = sprintf(
                                "The selected genes resolve to %s mapped features in %s. To avoid browser-sized plots, the expression plots are limited to the first %s mapped features. Download or inspect smaller gene batches when you need the omitted orthologs.",
                                format(resolution$feature_count_before_limit, big.mark = ","),
                                cross_integration_label(cross_key),
                                format(resolution$feature_limit, big.mark = ",")
                            ),
                            tone = "warning"
                        )
                    ))
                }

                if (length(genes) &&
                    (render_summary$panel_count > length(unique(genes)) ||
                        render_summary$mapped_feature_count > length(unique(genes)))) {
                    cards <- append(cards, list(
                        notice_card(
                            title = paste("Large orthogroup expansion in", cross_integration_label(cross_key)),
                            body = tagList(
                                tags$p(
                                    sprintf(
                                        "This view will render %s comparison blocks from %s selected source genes and %s mapped ortholog features. Large orthogroups can take longer to load, and very large mapped panels are trimmed by the rendering guardrail above.",
                                        format(render_summary$panel_count, big.mark = ","),
                                        format(length(unique(genes)), big.mark = ","),
                                        format(render_summary$mapped_feature_count, big.mark = ",")
                                    )
                                ),
                                if (length(render_summary$expansion_labels)) {
                                    tags$p(
                                        paste(
                                            "Largest expansions:",
                                            compact_gene_list(render_summary$expansion_labels, limit = 4),
                                            "."
                                        )
                                    )
                                }
                            ),
                            tone = "warning"
                        )
                    ))
                }

                div(class = "alert-stack", tagList(cards))
            })

            cross_mapping_table <- reactive({
                resolution <- cross_resolution()

                if (!nrow(resolution$plot_table)) {
                    return(tibble())
                }

                if (identical(integration_cfg$feature_mode, "medicago_space")) {
                    resolution$plot_table %>%
                        transmute(
                            `Source gene(s)` = source_gene_display,
                            `Medicago-space feature` = display_gene_labels(
                                "medicago",
                                target_feature_id,
                                include_gene_id_with_common = FALSE
                            ),
                            Orthogroup = ifelse(nzchar(orthogroup_label), orthogroup_label, "NA")
                        )
                } else {
                    resolution$plot_table %>%
                        transmute(
                            `Source gene(s)` = source_gene_display,
                            Species = target_species_label,
                            `SATURN feature` = target_display,
                            Orthogroup = ifelse(nzchar(orthogroup_label), orthogroup_label, "NA")
                        )
                }
            })

            output[[paste0(prefix, "_mapping_table_ui")]] <- renderUI({
                genes <- selected_source_genes()

                if (!length(genes)) {
                    return(div(
                        class = "summary-placeholder",
                        paste(
                            "Add source-species genes to see which features are used in the",
                            cross_integration_label(cross_key),
                            "tab."
                        )
                    ))
                }

                mapping_tbl <- cross_mapping_table()

                if (!nrow(mapping_tbl)) {
                    return(div(
                        class = "summary-placeholder",
                        paste(
                            "No selected genes resolve to the",
                            cross_integration_label(cross_key),
                            "feature set."
                        )
                    ))
                }

                html_summary_table(mapping_tbl)
            })

            ortholog_trace_source_genes <- reactive({
                genes <- selected_source_genes()

                if (!length(genes)) {
                    return(character(0))
                }

                genes
            })

            output[[paste0(prefix, "_ortholog_trace_notice_ui")]] <- renderUI({
                selected_n <- length(selected_source_genes())

                if (selected_n <= 10L) {
                    return(NULL)
                }

                div(
                    class = "alert-stack",
                    notice_card(
                        title = "Rendering the full ortholog trace",
                        body = sprintf(
                            "Showing all %d selected genes in the ortholog trace. Larger source panels and orthogroups now stay expanded here, so this section can take longer to render.",
                            selected_n
                        ),
                        tone = "warning"
                    )
                )
            })

            ortholog_trace_specs <- reactive({
                source_species <- current_comparison_source_species()
                source_genes <- ortholog_trace_source_genes()

                validate(
                    need(
                        length(source_genes) > 0,
                        "Add source-species genes to generate the ortholog trace panel."
                    )
                )

                lapply(source_genes, function(source_gene) {
                    source_label <- display_gene_labels(source_species, source_gene)[[1]]
                    orthogroup_id <- resolve_source_orthogroups(source_species, source_gene) %>%
                        filter(!is.na(orthogroup) & nzchar(orthogroup)) %>%
                        distinct(orthogroup) %>%
                        slice_head(n = 1) %>%
                        pull(orthogroup) %>%
                        first_nonempty()

                    species_specs <- lapply(within_species_keys, function(species_key) {
                        integration_method <- current_species_integration(species_key)

                        target_gene <- if (identical(species_key, source_species)) {
                            source_gene
                        } else if (is.na(orthogroup_id) || !nzchar(orthogroup_id)) {
                            NA_character_
                        } else {
                            first_nonempty(get_orthogroup_members(orthogroup_id, species_key))
                        }

                        target_label <- if (!is.na(target_gene) && nzchar(target_gene)) {
                            display_gene_labels(species_key, target_gene)[[1]]
                        } else {
                            NA_character_
                        }

                        if (is.na(target_gene) || !nzchar(target_gene)) {
                            return(list(
                                species_key = species_key,
                                integration_method = integration_method,
                                status = "missing_ortholog",
                                note = "No ortholog in this species",
                                target_gene = NA_character_,
                                target_label = NA_character_,
                                feature_id = NA_character_
                            ))
                        }

                        feature_id <- match_target_features(
                            target_species = species_key,
                            target_genes = target_gene,
                            integration_method = integration_method,
                            cross_space = FALSE
                        )[1] %||% NA_character_

                        if (is.na(feature_id) || !nzchar(feature_id)) {
                            return(list(
                                species_key = species_key,
                                integration_method = integration_method,
                                status = "missing_feature",
                                note = if (identical(species_key, source_species)) {
                                    "Gene not present in this atlas"
                                } else {
                                    "Ortholog not present in this atlas"
                                },
                                target_gene = target_gene,
                                target_label = target_label,
                                feature_id = NA_character_
                            ))
                        }

                        list(
                            species_key = species_key,
                            integration_method = integration_method,
                            status = "ok",
                            note = NA_character_,
                            target_gene = target_gene,
                            target_label = target_label,
                            feature_id = feature_id
                        )
                    })

                    list(
                        source_gene = source_gene,
                        source_label = source_label,
                        orthogroup = orthogroup_id,
                        title = if (!is.na(orthogroup_id) && nzchar(orthogroup_id)) {
                            paste0(source_label, " via orthogroup ", orthogroup_id)
                        } else {
                            paste0(source_label, " via orthogroup unavailable")
                        },
                        species_specs = species_specs
                    )
                })
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                ortholog_trace_source_genes(),
                current_species_integration("medicago"),
                current_species_integration("glycine"),
                current_species_integration("lotus"),
                cache = "app"
            )

            ortholog_trace_plot_obj <- reactive({
                row_specs <- ortholog_trace_specs()
                pt_size <- max(0.65, as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.45))

                row_plots <- lapply(row_specs, function(row_spec) {
                    expression_max <- max(unlist(lapply(row_spec$species_specs, function(spec) {
                        if (!identical(spec$status, "ok")) {
                            return(NA_real_)
                        }

                        obj <- get_within_object(spec$species_key, spec$integration_method)
                        expr_values <- aggregate_feature_expression(obj, spec$feature_id)

                        if (!length(expr_values)) {
                            return(NA_real_)
                        }

                        max(expr_values, na.rm = TRUE)
                    })), na.rm = TRUE)

                    if (!is.finite(expression_max) || expression_max <= 0) {
                        expression_max <- 1
                    }

                    species_panels <- lapply(row_spec$species_specs, function(spec) {
                        if (!identical(spec$status, "ok")) {
                            empty_panel <- empty_umap_message_plot(spec$note) +
                                theme(plot.margin = margin(4, 8, 6, 6))

                            return(add_species_caption(
                                empty_panel,
                                spec$species_key,
                                sublabel = spec$target_label
                            ))
                        }

                        obj <- get_within_object(spec$species_key, spec$integration_method)
                        feature_plot <- emphasized_feature_plot(
                            obj = obj,
                            feature_id = spec$feature_id,
                            fixed_max = expression_max,
                            pt_size = pt_size,
                            colorblind_safe = isTRUE(input$colorblind_safe),
                            expressing_size_boost = 0
                        ) +
                            labs(title = NULL, color = NULL) +
                            app_plot_theme() +
                            compact_feature_legend_guides() +
                            compact_feature_legend_theme() +
                            scale_x_continuous(expand = expansion(mult = 0.01)) +
                            scale_y_continuous(expand = expansion(mult = 0.01)) +
                            theme(
                                plot.title = element_blank(),
                                panel.grid = element_blank(),
                                axis.title = element_blank(),
                                axis.text = element_blank(),
                                axis.ticks = element_blank(),
                                axis.line = element_blank(),
                                plot.margin = margin(4, 8, 6, 6)
                            )

                        add_species_caption(
                            feature_plot,
                            spec$species_key,
                            sublabel = spec$target_label
                        )
                    })

                    wrap_titled_plot(
                        plot_obj = wrap_plots(
                            plotlist = species_panels,
                            ncol = 3,
                            guides = "collect"
                        ) & theme(legend.position = "top"),
                        title = row_spec$title
                    )
                })

                wrap_plots(plotlist = row_plots, ncol = 1)
            })

            cross_umap_plot_obj <- reactive({
                resolution <- cross_resolution()
                obj <- cross_object()
                panel_specs <- cross_comparison_panel_specs(resolution, cross_key)
                block_cols <- min(2L, max(1L, as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L)))
                pt_size <- max(0.65, as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.45))
                species_cells <- lapply(within_species_keys, function(species_key) {
                    cross_species_cell_ids(obj, species_key)
                })
                names(species_cells) <- within_species_keys

                validate(
                    need(
                        nrow(panel_specs) > 0,
                        paste(
                            "No selected genes resolve to features in the",
                            cross_integration_label(cross_key),
                            "integration."
                        )
                    )
                )

                reference_plot <- scCustomize::DimPlot_scCustom(
                    seurat_object = apply_metadata_display_order(obj, "species"),
                    colors_use = unname(distribution_color_map(obj@meta.data$species, "species")),
                    group.by = "species",
                    pt.size = max(0.45, pt_size * 0.85),
                    label = FALSE,
                    raster = TRUE
                ) &
                    app_plot_theme() &
                    theme(
                        legend.title = element_blank(),
                        legend.position = "none",
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        plot.margin = margin(4, 8, 6, 6)
                    )

                reference_panel <- wrap_titled_plot(
                    plot_obj = reference_plot + labs(title = NULL, color = NULL),
                    title = "Species overview"
                )

                comparison_blocks <- lapply(seq_len(nrow(panel_specs)), function(panel_idx) {
                    spec <- panel_specs[panel_idx, , drop = FALSE]
                    species_entries <- setNames(lapply(within_species_keys, function(species_key) {
                        spec[[paste0(species_key, "_panels")]][[1]]
                    }), within_species_keys)
                    species_expr <- setNames(lapply(within_species_keys, function(species_key) {
                        entries <- species_entries[[species_key]]

                        if (!panel_feature_count(entries)) {
                            return(list())
                        }

                        lapply(entries$feature_id, function(feature_id) {
                            aggregate_feature_expression(obj, feature_id)
                        })
                    }), within_species_keys)

                    expression_max <- max(unlist(lapply(within_species_keys, function(species_key) {
                        entries <- species_entries[[species_key]]
                        expr_list <- species_expr[[species_key]]
                        cell_ids <- species_cells[[species_key]]

                        if (!panel_feature_count(entries) || !length(cell_ids)) {
                            return(NA_real_)
                        }

                        vapply(seq_len(nrow(entries)), function(idx) {
                            expr_values <- expr_list[[idx]]

                            if (!length(expr_values)) {
                                return(NA_real_)
                            }

                            max(expr_values[cell_ids], na.rm = TRUE)
                        }, numeric(1))
                    })), na.rm = TRUE)

                    if (!is.finite(expression_max) || expression_max <= 0) {
                        expression_max <- 1
                    }

                    species_columns <- lapply(within_species_keys, function(species_key) {
                        entries <- species_entries[[species_key]]
                        expr_list <- species_expr[[species_key]]

                        if (!panel_feature_count(entries)) {
                            return(wrap_titled_plot(
                                plot_obj = empty_umap_message_plot("No mapped feature"),
                                title = species_label(species_key)
                            ))
                        }

                        feature_panels <- lapply(seq_len(nrow(entries)), function(idx) {
                            entry <- entries[idx, , drop = FALSE]

                            feature_plot <- emphasized_feature_plot(
                                obj = obj,
                                feature_id = entry$feature_id[[1]],
                                expression_values = expr_list[[idx]],
                                cell_ids = species_cells[[species_key]],
                                fixed_max = expression_max,
                                pt_size = pt_size,
                                colorblind_safe = isTRUE(input$colorblind_safe),
                                expressing_size_boost = 0
                            )

                            feature_plot <- feature_plot +
                                labs(title = NULL, color = NULL) +
                                app_plot_theme() +
                                compact_feature_legend_guides() +
                                compact_feature_legend_theme() +
                                scale_x_continuous(expand = expansion(mult = 0.01)) +
                                scale_y_continuous(expand = expansion(mult = 0.01)) +
                                theme(
                                    plot.title = element_blank(),
                                    panel.grid = element_blank(),
                                    axis.title = element_blank(),
                                    axis.text = element_blank(),
                                    axis.ticks = element_blank(),
                                    axis.line = element_blank(),
                                    plot.margin = margin(4, 8, 6, 6)
                                )

                            wrap_titled_plot(
                                plot_obj = feature_plot,
                                title = entry$feature_label[[1]],
                                size = 12,
                                hjust = 0.5
                            )
                        })

                        wrap_plots(plotlist = feature_panels, ncol = 1, guides = "collect") +
                            plot_title_annotation(species_label(species_key), size = 15, hjust = 0.5)
                    })

                    comparison_block <- wrap_plots(
                        plotlist = c(list(reference_panel), species_columns),
                        ncol = 4,
                        widths = c(0.72, 1, 1, 1),
                        guides = "collect"
                    ) &
                        theme(legend.position = "top")

                    comparison_block + plot_title_annotation(as.character(spec$title[[1]]))
                })

                wrap_plots(plotlist = comparison_blocks, ncol = block_cols)
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L),
                as.numeric(input[[paste0(prefix, "_pt_size")]] %||% 0.45),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            cross_dot_plot_obj <- reactive({
                resolution <- cross_resolution()
                obj <- cross_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste(
                            "No selected genes resolve to features in the",
                            cross_integration_label(cross_key),
                            "integration."
                        )
                    )
                )

                scCustomize::DotPlot_scCustom(
                    seurat_object = obj,
                    features = resolution$plot_features,
                    group.by = cross_group_by()
                ) +
                    scale_x_discrete(labels = resolution$label_map) +
                    app_plot_theme() +
                    labs(x = NULL, y = NULL, color = "Scaled expression", size = "% expressing") +
                    theme(
                        panel.grid.major.y = element_blank(),
                        axis.text.x = element_text(angle = 35, hjust = 1)
                    )
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                cross_group_by(),
                cache = "app"
            )

            cross_heatmap_plot_obj <- reactive({
                resolution <- cross_resolution()
                obj <- cross_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste(
                            "No selected genes resolve to features in the",
                            cross_integration_label(cross_key),
                            "integration."
                        )
                    )
                )

                build_expression_heatmap_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = cross_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                cross_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            cross_ridge_plot_obj <- reactive({
                resolution <- cross_resolution()
                obj <- cross_object()

                validate(
                    need(
                        length(resolution$plot_features) > 0,
                        paste(
                            "No selected genes resolve to features in the",
                            cross_integration_label(cross_key),
                            "integration."
                        )
                    )
                )

                build_expression_ridge_plot(
                    obj = obj,
                    feature_ids = resolution$plot_features,
                    label_map = resolution$label_map,
                    group_by = cross_group_by(),
                    colorblind_safe = isTRUE(input$colorblind_safe)
                )
            }) %>% bindCache(
                cross_key,
                input$source_species %||% "medicago",
                selected_source_genes(),
                cross_group_by(),
                isTRUE(input$colorblind_safe),
                cache = "app"
            )

            output[[paste0(prefix, "_umap_plot")]] <- renderPlot(
                {
                    cross_umap_plot_obj()
                },
                height = function() {
                    panel_specs <- tryCatch(cross_comparison_panel_specs(cross_resolution(), cross_key), error = function(e) tibble())
                    panel_n <- nrow(panel_specs)
                    block_cols <- min(2L, max(1L, as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L)))
                    row_units <- tryCatch({
                        block_units <- cross_comparison_height_units(panel_specs)

                        if (!length(block_units)) {
                            integer(0)
                        } else {
                            vapply(split(block_units, ceiling(seq_along(block_units) / block_cols)), max, integer(1))
                        }
                    }, error = function(e) integer(0))

                    if (!length(row_units)) {
                        return(560L)
                    }

                    max(560L, as.integer(180L + sum(row_units) * 310L))
                },
                res = 110
            )

            output[[paste0(prefix, "_ortholog_trace_plot")]] <- renderPlot(
                {
                    ortholog_trace_plot_obj()
                },
                height = function() {
                    ortholog_trace_height_px(length(ortholog_trace_source_genes()))
                },
                res = 110
            )

            output[[paste0(prefix, "_heatmap_plot")]] <- renderPlot(
                {
                    cross_heatmap_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
                    max(420, 110 + feature_n * 34)
                },
                res = 110
            )

            output[[paste0(prefix, "_ridge_plot")]] <- renderPlot(
                {
                    cross_ridge_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
                    group_n <- tryCatch(dplyr::n_distinct(cross_object()@meta.data[[cross_group_by()]]), error = function(e) 0L)
                    expression_ridge_height_px(feature_n = feature_n, group_n = group_n)
                },
                res = 110
            )

            output[[paste0(prefix, "_dot_plot")]] <- renderPlot(
                {
                    cross_dot_plot_obj()
                },
                height = function() {
                    feature_n <- tryCatch(length(cross_resolution()$plot_features), error = function(e) 0L)
                    max(420, 95 * feature_n + 120)
                },
                res = 110
            )

            output[[paste0("dl_", prefix, "_umap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_umap.", get_ext())
                },
                content = function(file) {
                    panel_specs <- cross_comparison_panel_specs(cross_resolution(), cross_key)
                    panel_n <- nrow(panel_specs)
                    block_cols <- min(2L, max(1L, as.integer(input[[paste0(prefix, "_umap_columns")]] %||% 1L)))
                    row_units <- cross_comparison_height_units(panel_specs)

                    if (length(row_units)) {
                        row_units <- vapply(split(row_units, ceiling(seq_along(row_units) / block_cols)), max, integer(1))
                    }

                    plot_height <- if (length(row_units)) {
                        max(6, 1.8 + sum(row_units) * 3.25)
                    } else {
                        6
                    }
                    plot_width <- if (block_cols == 1L) 17 else 30

                    save_ggplot(
                        file = file,
                        plot_obj = cross_umap_plot_obj(),
                        width = plot_width,
                        height = plot_height,
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "umap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_ortholog_trace")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_ortholog_trace.", get_ext())
                },
                content = function(file) {
                    gene_n <- length(ortholog_trace_source_genes())

                    save_ggplot(
                        file = file,
                        plot_obj = ortholog_trace_plot_obj(),
                        width = 14,
                        height = max(7, ortholog_trace_height_px(gene_n) / 95),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "ortholog_trace")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_heatmap")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_heatmap.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(cross_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = cross_heatmap_plot_obj(),
                        width = 13,
                        height = max(5.5, 1.8 + feature_n * 0.55),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "heatmap")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_ridge")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_ridgeplot.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(cross_resolution()$plot_features)
                    group_n <- dplyr::n_distinct(cross_object()@meta.data[[cross_group_by()]])

                    save_ggplot(
                        file = file,
                        plot_obj = cross_ridge_plot_obj(),
                        width = 10,
                        height = max(6, expression_ridge_height_px(feature_n, group_n) / 95),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "ridgeplot")
                    )
                }
            )

            output[[paste0("dl_", prefix, "_dot")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_dotplot.", get_ext())
                },
                content = function(file) {
                    feature_n <- length(cross_resolution()$plot_features)

                    save_ggplot(
                        file = file,
                        plot_obj = cross_dot_plot_obj(),
                        width = 12,
                        height = max(6, feature_n * 0.85 + 2),
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "dotplot")
                    )
                }
            )

            # ── Distribution UMAP ───────────────────────────────────────────
            cross_dist_umap_plot_obj <- reactive({
                obj <- cross_object()
                group_by <- cross_dist_group_by()
                pt_size <- as.numeric(input[[paste0(prefix, "_dist_pt_size")]] %||% 0.75)
                plot_obj <- apply_metadata_display_order(obj, group_by)
                colors_use <- distribution_colors_use(plot_obj, group_by)
                show_cluster_labels <- is_cluster_distribution_group(group_by)

                p <- scCustomize::DimPlot_scCustom(
                    seurat_object = plot_obj,
                    colors_use = colors_use,
                    group.by = group_by,
                    pt.size = pt_size,
                    label = show_cluster_labels,
                    repel = show_cluster_labels,
                    raster = TRUE
                ) &
                    app_plot_theme() &
                    theme(
                        legend.title = element_blank(),
                        legend.position = if (show_cluster_labels) "none" else "top",
                        panel.grid = element_blank(),
                        axis.title = element_blank(),
                        axis.text = element_blank(),
                        axis.ticks = element_blank(),
                        axis.line = element_blank(),
                        plot.margin = margin(8, 14, 10, 10)
                    )

                p & labs(title = NULL, color = NULL)
            }) %>% bindCache(
                cross_key,
                cross_dist_group_by(),
                as.numeric(input[[paste0(prefix, "_dist_pt_size")]] %||% 0.75),
                cache = "app"
            )

            output[[paste0(prefix, "_dist_umap_plot")]] <- renderPlot(
                { cross_dist_umap_plot_obj() },
                height = 540,
                res = 110
            )

            output[[paste0("dl_", prefix, "_dist_umap")]] <- downloadHandler(
                filename = function() paste0(prefix, "_dist_umap.", get_ext()),
                content = function(file) {
                    save_ggplot(
                        file = file,
                        plot_obj = cross_dist_umap_plot_obj(),
                        width = 10,
                        height = 8,
                        tab_label = cross_integration_label(cross_key),
                        integration_label = cross_integration_label(cross_key),
                        extra = list(atlas_export_type = "distribution_umap")
                    )
                }
            )

            cross_dist_umap3d_plot_data <- reactive({
                group_by <- cross_dist_group_by()
                obj <- apply_metadata_display_order(cross_object(), group_by)
                embedding <- get_cross_umap3d(cross_key)

                validate(need(!is.null(embedding), "3D UMAP is not available for this integration (no PCA or pre-computed 3D reduction found)."))

                cell_ids <- intersect(rownames(obj@meta.data), rownames(embedding))

                validate(need(length(cell_ids) > 0, "No 3D UMAP coordinates are available for this integration."))

                distribution_df <- tibble(
                    cell_id = cell_ids,
                    group_value = as.character(obj@meta.data[cell_ids, group_by, drop = TRUE]),
                    umap3d_1 = embedding[cell_ids, 1],
                    umap3d_2 = embedding[cell_ids, 2],
                    umap3d_3 = embedding[cell_ids, 3]
                ) %>%
                    filter(
                        !is.na(group_value) & nzchar(group_value),
                        !is.na(umap3d_1), !is.na(umap3d_2), !is.na(umap3d_3)
                    )

                validate(need(nrow(distribution_df) > 0, "No cells are available for the 3D UMAP view."))

                distribution_df <- distribution_df %>%
                    mutate(group_value = order_metadata_values(group_value, group_by))

                distribution_df <- stratified_point_sample(
                    distribution_df,
                    group_col = "group_value",
                    max_points = 30000L
                )

                full_color_map <- distribution_color_map(obj@meta.data[[group_by]], group_by)
                present_levels <- names(full_color_map)[names(full_color_map) %in% as.character(distribution_df$group_value)]

                list(
                    data = distribution_df,
                    color_map = full_color_map[present_levels],
                    group_by = group_by
                )
            })

            output[[paste0(prefix, "_dist_umap3d_plot")]] <- plotly::renderPlotly({
                plot_data <- cross_dist_umap3d_plot_data()
                df <- plot_data$data
                color_map <- plot_data$color_map
                group_by <- plot_data$group_by
                marker_size <- max(1.8, as.numeric(input[[paste0(prefix, "_dist_pt_size")]] %||% 0.75) * 2.2)
                group_levels <- names(color_map)

                p <- plotly::plot_ly(type = "scatter3d", mode = "markers")

                for (group_level in group_levels) {
                    group_df <- df %>% filter(as.character(group_value) == !!group_level)
                    if (!nrow(group_df)) next
                    p <- p %>%
                        plotly::add_trace(
                            data = group_df,
                            x = ~umap3d_1, y = ~umap3d_2, z = ~umap3d_3,
                            name = group_level,
                            marker = list(
                                size = marker_size,
                                color = unname(color_map[group_level]),
                                opacity = 0.8
                            ),
                            text = ~paste0(
                                metadata_column_label(group_by), ": ", group_value,
                                "<br>Cell: ", cell_id
                            ),
                            hovertemplate = "%{text}<extra></extra>"
                        )
                }

                p %>%
                    plotly::layout(
                        margin = list(l = 0, r = 180, b = 0, t = 10),
                        legend = list(
                            x = 1.02, xanchor = "left",
                            y = 1, yanchor = "top",
                            font = list(size = 14, color = unname(app_palette["text"])),
                            itemsizing = "constant",
                            bgcolor = "rgba(255,255,255,0.82)",
                            bordercolor = "rgba(201,214,196,0.92)",
                            borderwidth = 1
                        ),
                        scene = list(
                            aspectmode = "data",
                            xaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            yaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            zaxis = list(title = "", showgrid = FALSE, zeroline = FALSE, showticklabels = FALSE),
                            camera = list(eye = list(x = 1.22, y = 1.06, z = 0.75))
                        )
                    ) %>%
                    plotly::config(displaylogo = FALSE)
            })

            # ── Cluster composition ─────────────────────────────────────────
            cross_composition_by <- reactive({
                choices <- cross_composition_choices_cached()
                resolve_choice(
                    input[[paste0(prefix, "_composition_by")]],
                    choices,
                    default = if (length(choices)) unname(choices[[1]]) else NA_character_
                )
            })

            cross_composition_cluster_by <- reactive({
                clustering_cols <- c(
                    "Rank_1st_label", "Rank_2nd_label", "Rank_3rd_label", "Rank_4th_label", "Rank_5th_label",
                    "cluster_label", "species_cell_class",
                    "Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th"
                )
                available <- clustering_cols[clustering_cols %in% unname(cross_group_choices_cached())]
                preferred <- cross_dist_group_by()
                if (length(preferred) && preferred %in% available) return(preferred)
                if (length(available)) return(available[[1]])
                NA_character_
            })

            cross_markers_full <- reactive({
                read_cluster_markers_cache(cross_key, top_n = FALSE)
            })

            cross_marker_source <- reactive({
                marker_tbl <- cross_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NA_character_)
                }

                available_sources <- unique(as.character(marker_tbl$cluster_source))
                preferred_source <- cross_dist_group_by()
                rank_label_columns <- paste0(setdiff(distribution_cluster_columns, "cluster_label"), "_label")
                if (length(preferred_source) && preferred_source %in% rank_label_columns) {
                    preferred_source <- sub("_label$", "", preferred_source)
                }
                clustering_sources <- c("Rank_1st", "Rank_2nd", "Rank_3rd", "Rank_4th", "Rank_5th", "cluster_label", "species_cell_class")
                if (!(length(preferred_source) && preferred_source %in% clustering_sources)) {
                    preferred_source <- NULL
                }

                select_marker_cluster_source(
                    available_sources = available_sources,
                    preferred_source = preferred_source,
                    default_source = "Rank_1st"
                )
            })

            cross_marker_lookup <- reactive({
                get_cached_cluster_lookup(cross_key, cross_marker_source()) %||%
                    cluster_label_lookup(cross_object(), cross_marker_source())
            })

            cross_marker_choices <- reactive({
                marker_tbl <- cross_markers_full()

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(character(0))
                }

                marker_source <- cross_marker_source()

                marker_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source)

                available_clusters <- unique(as.character(marker_tbl$cluster))
                lookup_tbl <- cross_marker_lookup() %>%
                    filter(cluster %in% available_clusters)

                if (!nrow(lookup_tbl)) {
                    return(setNames(cluster_value_levels(available_clusters), cluster_value_levels(available_clusters)))
                }

                setNames(lookup_tbl$cluster, lookup_tbl$choice_label)
            })

            output[[paste0(prefix, "_marker_cluster_ui")]] <- renderUI({
                choices <- cross_marker_choices()

                if (!length(choices)) {
                    return(
                        tags$p(
                            class = "marker-status-hint",
                            tagList(
                                "Marker cache missing for this integration. Run ",
                                tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                                "."
                            )
                        )
                    )
                }

                selectInput(
                    inputId = paste0(prefix, "_marker_cluster"),
                    label = "Cluster",
                    choices = choices,
                    selected = resolve_choice(
                        input[[paste0(prefix, "_marker_cluster")]],
                        choices,
                        default = unname(choices[[1]])
                    )
                )
            })

            observeEvent(cross_marker_choices(), {
                pending_cluster <- pending_url_cluster()

                if (is.null(pending_cluster) || !nzchar(pending_cluster) || !identical(input$main_tabs %||% "overview", prefix)) {
                    return()
                }

                choices <- cross_marker_choices()

                if (!length(choices)) {
                    return()
                }

                selected_cluster <- resolve_choice(
                    pending_cluster,
                    choices,
                    default = unname(choices[[1]])
                )

                updateSelectInput(
                    session,
                    inputId = paste0(prefix, "_marker_cluster"),
                    selected = selected_cluster
                )
                pending_url_cluster(NULL)
            }, ignoreInit = FALSE)

            cross_marker_cluster <- reactive({
                choices <- cross_marker_choices()
                default_choice <- if (length(choices)) unname(choices[[1]]) else NULL

                resolve_choice(
                    input[[paste0(prefix, "_marker_cluster")]],
                    choices,
                    default = default_choice
                )
            })

            cross_marker_top_n <- reactive({
                marker_n <- suppressWarnings(as.integer(input[[paste0(prefix, "_marker_top_n")]] %||% 10L))

                if (is.na(marker_n) || marker_n < 1L) {
                    marker_n <- 10L
                }

                min(marker_n, 25L)
            })

            cross_marker_table_raw <- reactive({
                marker_tbl <- cross_markers_full()

                validate(
                    need(
                        !is.null(marker_tbl),
                        tagList(
                            "Cluster markers are not cached for this integration. Run ",
                            tags$code("Rscript scripts/build_cluster_markers_cache.R"),
                            "."
                        )
                    ),
                    need(nrow(marker_tbl) > 0, "No cluster markers are available for this integration.")
                )

                cluster_id <- cross_marker_cluster()
                marker_source <- cross_marker_source()
                filtered_tbl <- marker_tbl %>%
                    filter(cluster_source == !!marker_source) %>%
                    filter(cluster == !!cluster_id) %>%
                    arrange(p_val_adj, desc(avg_log2FC), desc(pct.1), gene)

                validate(need(nrow(filtered_tbl) > 0, "No marker rows are available for the selected cluster."))

                filtered_tbl
            })

            output[[paste0(prefix, "_markers_status_ui")]] <- renderUI({
                marker_tbl <- cross_markers_full()
                busy_message <- marker_job_messages()[[prefix]]

                if (is.null(marker_tbl) || !nrow(marker_tbl)) {
                    return(NULL)
                }

                if (!is.null(busy_message) && nzchar(busy_message)) {
                    return(tags$p(class = "marker-status-hint is-working", busy_message))
                }

                source_species <- current_comparison_source_species()
                cluster_choices <- cross_marker_choices()
                marker_source <- cross_marker_source()
                active_group_by <- cross_dist_group_by()
                current_cluster_label <- names(cluster_choices)[match(cross_marker_cluster(), unname(cluster_choices))]
                current_cluster_label <- current_cluster_label[!is.na(current_cluster_label) & nzchar(current_cluster_label)][1] %||% cross_marker_cluster()
                source_sentence <- if (marker_cluster_source_matches(marker_source, active_group_by)) {
                    paste0("using ", marker_cluster_source_label(marker_source), ".")
                } else {
                    paste0(
                        "using ",
                        marker_cluster_source_label(marker_source),
                        " because ",
                        metadata_column_label(active_group_by),
                        " is not a cached clustering solution."
                    )
                }
                tags$p(
                    class = "marker-status-hint",
                    sprintf(
                        "Showing markers for %s %s Top markers will be mapped back into the %s comparison panel and staged until you click Apply comparison panel.",
                        current_cluster_label,
                        source_sentence,
                        species_label(source_species)
                    )
                )
            })

            output[[paste0(prefix, "_markers_table")]] <- DT::renderDT({
                marker_tbl <- cross_marker_table_raw()

                display_tbl <- marker_tbl %>%
                    transmute(
                        gene = display_cross_feature_labels(cross_key, gene),
                        avg_log2FC = avg_log2FC,
                        pct.1 = pct.1,
                        pct.2 = pct.2,
                        p_val_adj = p_val_adj
                    )

                DT::datatable(
                    display_tbl,
                    rownames = FALSE,
                    escape = TRUE,
                    selection = "none",
                    class = "stripe hover order-column compact",
                    options = list(
                        dom = "tip",
                        pageLength = 10,
                        autoWidth = TRUE,
                        order = list(list(1, "desc"))
                    )
                ) %>%
                    DT::formatRound(columns = c("avg_log2FC", "pct.1", "pct.2"), digits = 3) %>%
                    DT::formatSignif(columns = "p_val_adj", digits = 3)
            })

            observeEvent(input[[paste0(prefix, "_add_markers")]], {
                button_id <- paste0(prefix, "_add_markers")
                marker_tbl <- cross_marker_table_raw() %>%
                    slice_head(n = cross_marker_top_n())

                source_species <- current_comparison_source_species()
                set_marker_job_message(
                    prefix,
                    sprintf(
                        "Adding the top %d markers from the %s cluster in %s into the shared %s comparison panel...",
                        nrow(marker_tbl),
                        cross_marker_cluster(),
                        cross_integration_label(cross_key),
                        species_label(source_species)
                    )
                )
                session$sendCustomMessage("atlas_button_busy", list(
                    button_id = button_id,
                    label = "Adding..."
                ))
                on.exit(
                    session$onFlushed(function() {
                        session$sendCustomMessage("atlas_button_busy", list(
                            button_id = button_id,
                            busy = FALSE
                        ))
                    }, once = TRUE),
                    add = TRUE
                )
                candidate_genes <- map_cross_marker_features_to_source_genes(
                    cross_key = cross_key,
                    feature_ids = marker_tbl$gene,
                    source_species = source_species
                )

                valid_source_ids <- build_gene_choices(
                    source_species,
                    source_integration()
                )$feature_ids %||% character(0)

                genes_to_add <- unique(candidate_genes[candidate_genes %in% valid_source_ids])

                if (!length(genes_to_add)) {
                    showNotification(
                        sprintf(
                            "No top markers from this %s cluster could be mapped into the %s comparison panel.",
                            cross_integration_label(cross_key),
                            species_label(source_species)
                        ),
                        type = "warning",
                        duration = 8
                    )
                    return()
                }

                updated_selection <- unique(c(staged_source_genes(), genes_to_add))

                update_selected_genes_input(
                    choice_bundle = build_gene_choices(
                        source_species,
                        source_integration()
                    ),
                    selected = updated_selection
                )

                skipped_n <- nrow(marker_tbl) - length(genes_to_add)
                showNotification(
                    paste0(
                        "Added ",
                        length(genes_to_add),
                        " mapped marker gene(s) to the shared selection. Click Apply comparison panel to refresh the cross-species views",
                        if (skipped_n > 0) paste0("; ", skipped_n, " could not be mapped into the current source species.") else "."
                    ),
                    type = "message",
                    duration = 6
                )
            }, ignoreInit = TRUE)

            output[[paste0("dl_", prefix, "_markers")]] <- downloadHandler(
                filename = function() {
                    paste0(prefix, "_cluster_", cross_marker_cluster(), "_markers.csv")
                },
                content = function(file) {
                    marker_tbl <- cross_marker_table_raw()
                    lookup_tbl <- cross_marker_lookup()
                    cluster_labels <- lookup_tbl$cluster_label[match(marker_tbl$cluster, lookup_tbl$cluster)]

                    export_tbl <- marker_tbl %>%
                        mutate(
                            cluster_label = ifelse(is.na(cluster_labels) | !nzchar(cluster_labels), cluster, cluster_labels),
                            gene_label = display_cross_feature_labels(cross_key, gene)
                        ) %>%
                        select(cluster, cluster_label, gene, gene_label, avg_log2FC, pct.1, pct.2, p_val_adj) %>%
                        add_export_provenance_columns(
                            tab_label = cross_integration_label(cross_key),
                            integration_label = cross_integration_label(cross_key),
                            extra = list(atlas_export_type = "cluster_markers")
                        )

                    write.csv(export_tbl, file = file, row.names = FALSE, na = "")
                }
            )

            cross_composition_plot_obj <- reactive({
                obj <- cross_object()
                composition_by <- cross_composition_by()
                cluster_by <- cross_composition_cluster_by()
                plot_obj <- apply_metadata_display_order(obj, composition_by)

                validate(
                    need(!is.null(composition_by) && nzchar(composition_by %||% ""), "No composition metadata available for this integration."),
                    need(!is.na(cluster_by) && nzchar(cluster_by), "No clustering metadata available for this integration.")
                )

                md <- plot_obj@meta.data
                composition_df <- tibble(
                    cluster = as.character(md[[cluster_by]]),
                    composition = as.character(md[[composition_by]])
                ) %>%
                    filter(
                        !is.na(cluster) & nzchar(cluster),
                        !is.na(composition) & nzchar(composition)
                    ) %>%
                    count(cluster, composition, name = "cell_count")

                validate(need(nrow(composition_df) > 0, "No cluster composition data available."))

                composition_df <- composition_df %>%
                    mutate(
                        cluster = factor(cluster, levels = cluster_value_levels(cluster)),
                        composition = order_metadata_values(composition, composition_by)
                    )

                fill_values <- composition_colors_use(composition_df$composition, composition_by)
                legend_rows <- max(1L, min(3L, ceiling(dplyr::n_distinct(composition_df$composition) / 8L)))

                ggplot(composition_df, aes(x = cluster, y = cell_count, fill = composition)) +
                    geom_col(position = "fill", width = 0.82, colour = "#FFFFFF", linewidth = 0.18) +
                    scale_y_continuous(
                        labels = scales::label_percent(accuracy = 1),
                        expand = expansion(mult = c(0, 0.02))
                    ) +
                    {if (!is.null(fill_values)) scale_fill_manual(values = fill_values)} +
                    guides(fill = guide_legend(nrow = legend_rows, byrow = TRUE)) +
                    labs(x = metadata_column_label(cluster_by), y = "Percent of cells", fill = NULL) +
                    app_plot_theme() +
                    theme(
                        panel.grid.major.x = element_blank(),
                        axis.text.x = element_text(size = 11, angle = 35, hjust = 1),
                        legend.position = "top",
                        legend.text = element_text(size = 9),
                        plot.margin = margin(8, 12, 8, 10)
                    )
            }) %>% bindCache(
                cross_key,
                cross_composition_by(),
                cross_composition_cluster_by(),
                cache = "app"
            )

            output[[paste0(prefix, "_composition_plot")]] <- renderPlot(
                { cross_composition_plot_obj() },
                height = function() {
                    composition_by <- tryCatch(cross_composition_by(), error = function(e) NULL)
                    obj <- tryCatch(cross_object(), error = function(e) NULL)
                    if (is.null(obj) || is.null(composition_by) || !nzchar(composition_by %||% "")) return(460)
                    level_count <- tryCatch({
                        vals <- as.character(obj@meta.data[[composition_by]])
                        length(unique(vals[!is.na(vals) & nzchar(vals)]))
                    }, error = function(e) 1L)
                    legend_rows <- max(1L, min(3L, ceiling(level_count / 8L)))
                    max(460, 380 + legend_rows * 44)
                },
                res = 110
            )
        })
    }

    walk(cross_integration_keys, register_cross_tab)
}

shinyApp(ui = ui, server = server)
