namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ設定
    /// </summary>
    internal class GraphConfig
    {
        /// <summary>
        /// グラフ設定を作成する。
        /// </summary>
        /// <param name="xConfig">X軸設定</param>
        /// <param name="yConfig">Y軸設定</param>
        /// <param name="option">グラフサイズOption</param>
        public GraphConfig(XAxisConfig xConfig, YAxisConfig yConfig, GraphOption option)
        {
            XAxisConfig = xConfig;
            YAxisConfig = yConfig;
            FigureWidth = option.FigureWidth;
            FigureHeight = option.FigureHeight;
        }

        /// <summary>
        /// X軸設定を取得する。
        /// </summary>
        public XAxisConfig XAxisConfig { get; }

        /// <summary>
        /// Y軸設定を取得する。
        /// </summary>
        public YAxisConfig YAxisConfig { get; }

        /// <summary>
        /// グラフ幅を取得する。
        /// </summary>
        public FigureWidth FigureWidth { get; }

        /// <summary>
        /// グラフの高さを取得する。
        /// </summary>
        public FigureHeight FigureHeight { get; }
    }
}
