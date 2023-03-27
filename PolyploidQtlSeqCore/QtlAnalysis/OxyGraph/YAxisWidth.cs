namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// Y軸幅
    /// </summary>
    internal class YAxisWidth
    {
        private readonly YAxisConfig _yAxisConfig;

        /// <summary>
        /// Y軸幅インスタンスを作成する。
        /// </summary>
        /// <param name="yAxis">Y軸</param>
        /// <param name="width">Y軸幅(pixel)</param>
        public YAxisWidth(YAxisConfig yAxis, double width)
        {
            _yAxisConfig = yAxis;
            Width = width;
        }

        /// <summary>
        /// Y軸幅を取得する。
        /// </summary>
        public double Width { get; }

        /// <summary>
        /// Y軸設定調整を行う。
        /// </summary>
        /// <param name="maxWidth">最大幅</param>
        /// <returns>調整済みYAxisConfig</returns>
        public YAxisConfig Correct(YAxisWidth maxWidth)
        {
            if (this == maxWidth) return _yAxisConfig;
            if (Width == maxWidth.Width) return _yAxisConfig;

            var diff = maxWidth.Width - Width;
            return _yAxisConfig.Correct(diff);
        }
    }
}
