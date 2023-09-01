namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ設定
    /// </summary>
    internal class GraphSettings
    {
        private static readonly IReadOnlyDictionary<string, string> _toLongNameDictionary = new Dictionary<string, string>()
        {
            [FigureWidth.SHORT_NAME] = FigureWidth.LONG_NAME,
            [FigureWidth.LONG_NAME] = FigureWidth.LONG_NAME,

            [FigureHeight.SHORT_NAME] = FigureHeight.LONG_NAME,
            [FigureHeight.LONG_NAME] = FigureHeight.LONG_NAME,

            [XAxisMajorStep.SHORT_NAME] = XAxisMajorStep.LONG_NAME,
            [XAxisMajorStep.LONG_NAME] = XAxisMajorStep.LONG_NAME
        };


        /// <summary>
        /// LongName変換項目を辞書に追加する。
        /// </summary>
        /// <param name="dictionary">LongName変換辞書</param>
        public static void AddLongNameKeyValuePair(Dictionary<string, string> dictionary)
        {
            foreach (var keyValuePair in _toLongNameDictionary)
            {
                dictionary.Add(keyValuePair.Key, keyValuePair.Value);
            }
        }

        /// <summary>
        /// グラフ設定を作成する。
        /// </summary>
        /// <param name="settingValue">設定値値</param>
        /// <param name="parameterDictionary">LongNameパラメーター辞書</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public GraphSettings(IGraphSettingValue settingValue, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            FigureWidth = new FigureWidth(settingValue.FigureWidth, parameterDictionary, userOptionDictionary);
            FigureHeight = new FigureHeight(settingValue.FigureHeight, parameterDictionary, userOptionDictionary);
            XAxisMajorStep = new XAxisMajorStep(settingValue.XAxisMajorStep, parameterDictionary, userOptionDictionary);
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

        /// <summary>
        /// パラメータファイルに記載する行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        public string[] ToParameterFileLines()
        {
            return new[]
            {
                FigureWidth.ToParameterFileLine(),
                FigureHeight.ToParameterFileLine(),
                XAxisMajorStep.ToParameterFileLine()
            };
        }
    }
}
