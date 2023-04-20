using PolyploidQtlSeqCore.Options;

namespace PolyploidQtlSeqCore.QtlAnalysis.Distribution
{
    /// <summary>
    /// 親2のPlex数
    /// </summary>
    public class Parent2PlexNumber
    {
        /// <summary>
        /// plex数の最小値
        /// </summary>
        public const int MINIMUM = 1;

        /// <summary>
        /// plex数の規定値
        /// </summary>
        public const int DEFAULT = 1;

        /// <summary>
        /// plex数の仮の最大値
        /// </summary>
        public const int MAXIMUM = 100;

        /// <summary>
        /// オプションスイッチのShortName
        /// </summary>
        public const string SHORT_NAME = "np";

        /// <summary>
        /// オプションスイッチのLongName
        /// </summary>
        public const string LONG_NAME = "NPlex";

        /// <summary>
        /// オプションスイッチの説明
        /// </summary>
        public const string DESCRIPTION = "Specify the plexity of Parent2 used for QTL analysis.";

        /// <summary>
        /// データ検証エラーメッセージ
        /// </summary>
        public const string VALIDATION_ERROR_MESSAGE = "The -np option must be an integer greater than or equal to 1 and less than or equal to ploidy.";

        /// <summary>
        /// 親2のPlex数を作成する。
        /// </summary>
        /// <param name="plex">plex数</param>
        /// <param name="parameterDictionary">パラメータファイルの中身</param>
        /// <param name="userOptionDictionary">ユーザー指定LongName辞書</param>
        public Parent2PlexNumber(int plex, IReadOnlyDictionary<string, string> parameterDictionary,
            IReadOnlyDictionary<string, bool> userOptionDictionary)
        {
            Value = OptionValue.GetValue(LONG_NAME, plex, parameterDictionary, userOptionDictionary);

            // 最大値はploidyによって変わるのでここでは判断できない
            if (Value < MINIMUM || Value > MAXIMUM) throw new ArgumentException(VALIDATION_ERROR_MESSAGE);
        }

        /// <summary>
        /// Plex数を取得する。
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
