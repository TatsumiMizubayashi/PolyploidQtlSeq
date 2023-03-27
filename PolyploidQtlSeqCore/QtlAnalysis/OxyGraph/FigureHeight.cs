using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ画像の高さ
    /// </summary>
    public class FigureHeight
    {
        /// <summary>
        /// グラフ高さの最小値
        /// </summary>
        public const int MINIMUM = 100;

        /// <summary>
        /// グラフ高さの最大値
        /// </summary>
        public const int MAXIMUM = 2000;

        /// <summary>
        /// グラフ高さの規定値
        /// </summary>
        public const int DEFAULT = 300;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "fh";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "figHeight";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Graph image height (Pixel).";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -fh option must be an integer greater than 100 and less than or equal to 2000.";

        /// <summary>
        /// グラフ画像の高さ(pixel)を作成する。
        /// </summary>
        /// <param name="height">高さ(pixel)</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public FigureHeight(int height, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, height, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// グラフ画像の高さ(Pixel)を取得する。
        /// </summary>
        internal int Value { get; }

        /// <summary>
        /// パラメータファイル記載用行テキストに変換する。
        /// </summary>
        /// <returns>パラメータ行テキスト</returns>
        internal string ToParameterFileLine()
        {
            return $"{LONG_NAME}\t{Value}";
        }
    }
}
