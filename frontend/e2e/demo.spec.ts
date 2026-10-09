import { test, expect } from '@playwright/test';

test('actual synthetic geometry, charts, linked selection and camera interaction', async ({ page, request }) => {
  const errors: string[] = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto('/');
  await expect(page.getByLabel('Case', { exact: true })).toHaveValue('synthetic_demo');
  await expect(page.getByRole('note')).toContainText('stylized generated tubes');
  const viewer = page.locator('.vtk-wrapper');
  await expect(viewer).toHaveAttribute('data-loaded-count', '12');
  await expect(page.locator('.vessel-button')).toHaveCount(12);
  await expect(page.locator('.recharts-line-curve')).toHaveCount(2);
  const profile = await (await request.get('http://127.0.0.1:8000/api/cases/synthetic_demo/profiles/1')).json();
  await page.getByRole('button', { name: 'Show minimum caliber' }).click();
  const sample = profile.samples.reduce((a: {distance_mm:number}, b: {distance_mm:number}) => Math.abs(b.distance_mm-profile.minimum_location_mm)<Math.abs(a.distance_mm-profile.minimum_location_mm) ? b : a);
  await expect(viewer).toHaveAttribute('data-marker', JSON.stringify([sample.x_mm,sample.y_mm,sample.z_mm]));
  await expect(page.locator('.profile-selection')).toContainText('Selected location');
  await page.getByRole('button', { name: 'SYN-02 Synthetic tube 02' }).click();
  await expect(viewer).toHaveAttribute('data-marker','none');
  await expect(page.locator('.profile-chart-title')).toContainText('Synthetic tube 02');
  const chart = page.locator('.recharts-wrapper');
  await chart.scrollIntoViewIfNeeded();
  const box = (await chart.boundingBox())!;
  await page.mouse.move(box.x+box.width*0.55,box.y+box.height*0.5);
  await page.mouse.click(box.x+box.width*0.55,box.y+box.height*0.5);
  await expect(viewer).not.toHaveAttribute('data-marker','none');
  const canvas = page.locator('.vtk-container canvas');
  const before = await canvas.screenshot();
  const initialCamera = await page.locator(".vtk-container").getAttribute("data-camera");
  const area = (await canvas.boundingBox())!;
  await page.mouse.move(area.x+area.width/2,area.y+area.height/2);
  await page.mouse.down();
  await page.mouse.move(area.x+area.width/2+70,area.y+area.height/2+40,{steps:10});
  await page.mouse.up();
  expect((await canvas.screenshot()).equals(before)).toBeFalsy();
  const rotated = await canvas.screenshot();
  await page.mouse.wheel(0,-200);
  await page.waitForTimeout(300);
  expect((await canvas.screenshot()).equals(rotated)).toBeFalsy();
  await page.getByRole('button',{name:'Reset view'}).click();
  await page.waitForTimeout(300);
  expect(JSON.parse((await page.locator(".vtk-container").getAttribute("data-camera"))!)).toEqual(JSON.parse(initialCamera!));
  expect(errors).toEqual([]);
});

test('prepared legacy case switching clears selection and explicit missing cases fail', async ({page,request}) => {
  const report = await (await request.get('http://127.0.0.1:8000/api/cases/synthetic_demo')).json();
  // A legacy fixture is transported through API mocks; actual synthetic assets remain real.
  await page.route('**/api/cases', route => route.fulfill({json:[{case_id:'synthetic_demo'},{case_id:'legacy_case'}]}));
  await page.route('**/api/cases/legacy_case', route => route.fulfill({json:{...report,case_id:'legacy_case',provenance:undefined}}));
  await page.route('**/api/cases/legacy_case/geometry/*', async route => route.fulfill({response:await request.get(route.request().url().replace('legacy_case','synthetic_demo'))}));
  await page.route('**/api/cases/legacy_case/profiles/*', async route => {
    const profile = await (await request.get(route.request().url().replace('legacy_case','synthetic_demo'))).json();
    await route.fulfill({json:{...profile,case_id:'legacy_case'}});
  });
  await page.goto('/');
  await page.getByRole('button',{name:'Show minimum caliber'}).click();
  await page.getByLabel('Case',{exact:true}).selectOption('legacy_case');
  await expect(page.locator('.vtk-wrapper')).toHaveAttribute('data-marker','none');
  await expect(page.locator('.vtk-wrapper')).toHaveAttribute('data-loaded-count','12');
  await expect(page.getByRole('note')).toHaveCount(0);
  await expect(page).toHaveURL(/case=legacy_case/);
  await page.goto('/?case=missing_real_case');
  await expect(page.locator('.status-screen')).toContainText('Unable to load case missing_real_case');
  await expect(page.locator('.metrics')).toHaveCount(0);
});

test('empty cases and empty profiles are explicit',async ({page,request}) => {
  await page.route('**/api/cases',route=>route.fulfill({json:[]}));
  await page.goto('/');
  await expect(page.locator('.status-screen')).toContainText('No prepared cases');
  const report = await (await request.get('http://127.0.0.1:8000/api/cases/synthetic_demo')).json();
  await page.route('**/api/cases/empty',route=>route.fulfill({json:{...report,vessels:[],vessel_count:0}}));
  await page.goto('/?case=empty');
  await expect(page.locator('.status-screen')).toContainText('Case has no vessel metrics');
  await page.route('**/profiles/*',route=>route.fulfill({json:{case_id:'synthetic_demo',label:1,name:'Synthetic tube 01',samples:[]}}));
  await page.goto('/?case=synthetic_demo');
  await expect(page.locator('.profile-status')).toContainText('No caliber samples available');
});
