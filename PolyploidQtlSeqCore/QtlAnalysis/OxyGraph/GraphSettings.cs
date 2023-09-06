namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ設定
    /// </summary>
    internal class GraphSettings
    {
        /// <summary>
        /// グラフ設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public GraphSettings(IGraphSettingValue settingValue)
        {
            FigureWidth = new FigureWidth(settingValue.FigureWidth);
            FigureHeight = new FigureHeight(settingValue.FigureHeight);
            XAxisMajorStep = new XAxisMajorStep(settingValue.XAxisMajorStep);
        }

        /// <summary>
        /// グラフ画像の幅を取得する。
        /// </summary>
        public FigureWidth FigureWidth { get; }

        /// <summary>
        /// グラフ画像の高さを取得する。
        /// </summary>
        public FigureHeight FigureHeight { get; }

        /// <summary>
        /// X軸目盛り間隔(MB)を取得する。
        /// </summary>
        public XAxisMajorStep XAxisMajorStep { get; }

    }
}
