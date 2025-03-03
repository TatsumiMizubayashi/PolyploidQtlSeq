using OxyPlot.Series;
using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフデータ
    /// </summary>
    internal class GraphData
    {
        private readonly AllVariants _allVariants;
        private readonly AllWindows _allWindows;

        /// <summary>
        /// グラフデータを作成する。
        /// </summary>
        /// <param name="allVariants">全変異情報</param>
        /// <param name="allWindows">全Window</param>
        /// <param name="majorStep">X軸MajorStep</param>
        public GraphData(SnpIndexVariantWithSlidingWindowQtl[] allVariants, Window[] allWindows, XAxisMajorStep majorStep)
        {
            _allVariants = new AllVariants(allVariants);
            _allWindows = new AllWindows(allWindows);

            ChrNames = _allVariants.GetChrNames();

            var xAxis = _allWindows.CreateXAsisConfig(majorStep);
            var snpIndexYAxis = _allVariants.CreateBulkSnpIndexYAxisConfig();
            var deltaSnpIndexYAxis = DeltaSnpIndexYAxisConfigCreator.Create(allVariants, allWindows);
            var scoreYAxis = _allWindows.CreateWindowScoreYAxisConfig();
            var qtlCountYAsis = _allWindows.CreateQtlCountYAxisConfig();
            GraphAxes = new GraphAxes(xAxis, snpIndexYAxis, deltaSnpIndexYAxis, scoreYAxis, qtlCountYAsis);
        }

        /// <summary>
        /// 染色体名を取得する。
        /// </summary>
        public string[] ChrNames { get; }

        /// <summary>
        /// グラフ軸を取得する。
        /// </summary>
        public GraphAxes GraphAxes { get; }

        /// <summary>
        /// 指定した染色体のBulk1 SnpIndex ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>Bulk1 SnpIndex ScatterSeries</returns>
        public ScatterSeries CreateBulk1SnpIndexScatterSeries(string chrName)
        {
            var chrVariant = _allVariants.GetVariants(chrName);
            return chrVariant.CreateBulk1SnpIndexScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のBulk1 SnpIndex=0 ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>Bulk1 SnpIndex=0 ScatterSeries</returns>
        public ScatterSeries CreateBulk1ZeroSnpIndexScatterSeries(string chrName)
        {
            var chrVariant = _allVariants.GetVariants(chrName);
            return chrVariant.CreateBulk1ZeroSnpIndexScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のBulk2 SnpIndex ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>Bulk2 SnpIndex ScatterSeries</returns>
        public ScatterSeries CreateBulk2SnpIndexScatterSeries(string chrName)
        {
            var chrVariant = _allVariants.GetVariants(chrName);
            return chrVariant.CreateBulk2SnpIndexScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のBulk2 SnpIndex=0 ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>Bulk2 SnpIndex=0 ScatterSeries</returns>
        public ScatterSeries CreateBulk2ZeroSnpIndexScatterSeries(string chrName)
        {
            var chrVariant = _allVariants.GetVariants(chrName);
            return chrVariant.CreateBulk2ZeroSnpIndexScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のΔSnpIndex ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>ΔSnpIndex ScatterSeries</returns>
        public ScatterSeries CreateDeltaSnpIndexScatterSeries(string chrName)
        {
            var chrVariant = _allVariants.GetVariants(chrName);
            return chrVariant.CreateDeltaSnpIndexScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のスコア ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>スコア ScatterSeries</returns>
        public ScatterSeries CreateScoreScatterSeries(string chrName)
        {
            var chrWindows = _allWindows.GetZeroExcludeWindows(chrName);
            return chrWindows.CreateScorexScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のP95QTL数 ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>P95QTL数 ScatterSeries</returns>
        public ScatterSeries CreateP95QtlCountScatterSeries(string chrName)
        {
            var chrWindows = _allWindows.GetZeroExcludeWindows(chrName);
            return chrWindows.CreateP95QtlCountScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のP99QTL数 ScatterSeriesを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>P99QTL数 ScatterSeries</returns>
        public ScatterSeries CreateP99QtlCountScatterSeries(string chrName)
        {
            var chrWindows = _allWindows.GetZeroExcludeWindows(chrName);
            return chrWindows.CreateP99QtlCountScatterSeries();
        }

        /// <summary>
        /// 指定した染色体のBulk1平均SnpIndex LineSeriesを取得する。
        /// 個体数0Windowがある場合は分割され新しいLineSeriesになる。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>Bulk1平均SnpIndex LineSeries</returns>
        public LineSeries[] CreateBulk1AverageSnpIndexLineSeries(string chrName)
        {
            var chrWindowsList = _allWindows.GetZeroSplitWindowsList(chrName);
            return [.. chrWindowsList.Select(x => x.CreateBulk1AverageSnpIndexLineSeries())];
        }

        /// <summary>
        /// 指定した染色体のBulk2平均SnpIndex LineSeriesを取得する。
        /// 個体数0Windowがある場合は分割され新しいLineSeriesになる。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>Bulk2平均SnpIndex LineSeries</returns>
        public LineSeries[] CreateBulk2AverageSnpIndexLineSeries(string chrName)
        {
            var chrWindowsList = _allWindows.GetZeroSplitWindowsList(chrName);
            return [.. chrWindowsList.Select(x => x.CreateBulk2AverageSnpIndexLineSeries())];
        }


        /// <summary>
        /// 指定した染色体の平均ΔSnpIndex LineSeriesを取得する。
        /// 個体数0Windowがある場合は分割され新しいLineSeriesになる。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>平均ΔSnpIndex LineSeries</returns>
        public LineSeries[] CreateAverageDeltaSnpIndexLineSeries(string chrName)
        {
            var chrWindowsList = _allWindows.GetZeroSplitWindowsList(chrName);
            return [.. chrWindowsList.Select(x => x.CreateAverageDeltaSnpIndexLineSeries())];
        }

        /// <summary>
        /// 指定した染色体のP95しきい値 LineSeriesを取得する。
        /// 個体数0Windowがある場合は分割され新しいLineSeriesになる。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>P95しきい値 LineSeries</returns>
        public LineSeries[] CreateP95ThresholdLineSeries(string chrName)
        {
            var chrWindowsList = _allWindows.GetZeroSplitWindowsList(chrName);
            return [.. chrWindowsList
                .SelectMany(x =>
                {
                    var (plus, minus) = x.CreateP95ThresholdLineSeries();
                    return new[] { plus, minus };
                })];
        }

        /// <summary>
        /// 指定した染色体のP99しきい値 LineSeriesを取得する。
        /// 個体数0Windowがある場合は分割され新しいLineSeriesになる。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>P99しきい値 LineSeries</returns>
        public LineSeries[] CreateP99ThresholdLineSeries(string chrName)
        {
            var chrWindowsList = _allWindows.GetZeroSplitWindowsList(chrName);
            return [.. chrWindowsList
                .SelectMany(x =>
                {
                    var (plus, minus) = x.CreateP99ThresholdLineSeries();
                    return new[] { plus, minus };
                })];
        }

    }
}
