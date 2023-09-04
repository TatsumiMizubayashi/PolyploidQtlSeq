using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// グラフ画像の幅
    /// </summary>
    public class FigureWidth
    {
        /// <summary>
        /// グラフ幅の最小値
        /// </summary>
        public const int MINIMUM = 300;

        /// <summary>
        /// グラフ幅の最大値
        /// </summary>
        public const int MAXIMUM = 5000;

        /// <summary>
        /// グラフ幅の規定値
        /// </summary>
        [Obsolete("削除予定")]
        public const int DEFAULT = 1200;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        [Obsolete("削除予定")]
        public const string SHORT_NAME = "fw";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        [Obsolete("削除予定")]
        public const string LONG_NAME = "figWidth";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        [Obsolete("削除予定")]
        public const string DESCRIPTION = "Width (pixel) of the graph images.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        [Obsolete("削除予定")]
        public const string VALIDATION_ERROR_MESSAGE = "The -fw option must be an integer greater than 300 and less than or equal to 5000.";

        /// <summary>
        /// グラフ画像の幅を作成する。
        /// </summary>
        /// <param name="width">幅(pixel)</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public FigureWidth(int width, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, width, parameterDictionary, userOptionDictionary);

            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// グラフ画像の幅(Pixel)を取得する。
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
